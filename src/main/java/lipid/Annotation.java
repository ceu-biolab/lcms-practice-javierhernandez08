package lipid;

import adduct.Adduct;
import adduct.AdductList;

import java.util.Collections;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IoniationMode ionizationMode;
    private String adduct; // The adduct will be detected based on the groupedSignals
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;

    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode, Collections.emptySet());
    }

    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
        this.groupedSignals = new TreeSet<>(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;

        if (!this.groupedSignals.isEmpty()) {
            detectAdductFromGroupedSignals();
        }
    }

    private void detectAdductFromGroupedSignals() {
        double toleranceDa = 0.01;
        double waterLossMass = 18.0106;
        Peak[] peaks = groupedSignals.toArray(new Peak[0]);

        for (int i = 0; i < peaks.length; i++) {
            double mz1 = peaks[i].getMz();

            for (String adduct1 : AdductList.MAPMZPOSITIVEADDUCTS.keySet()) {
                Double mass1 = Adduct.getMonoisotopicMassFromMZ(mz1, adduct1);
                if (mass1 == null) continue;

                for (int j = i + 1; j < peaks.length; j++) {
                    double mz2 = peaks[j].getMz();

                    for (String adduct2 : AdductList.MAPMZPOSITIVEADDUCTS.keySet()) {
                        if (adduct1.equals(adduct2)) continue;

                        Double mass2 = Adduct.getMonoisotopicMassFromMZ(mz2, adduct2);
                        if (mass2 == null) continue;

                        if (Math.abs(mass1 - mass2) <= toleranceDa) {
                            this.setAdduct(mz1 < mz2 ? adduct1 : adduct2);
                            return;
                        }
                    }

                    if (Math.abs(mz1 - mz2 - waterLossMass) <= toleranceDa || Math.abs(mz2 - mz1 - waterLossMass) <= toleranceDa) {
                        double mzMayor = Math.max(mz1, mz2);
                        String adductMayor = findBestMatchingAdduct(mzMayor);
                        if (adductMayor != null) {
                            this.setAdduct(adductMayor);
                            return;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < peaks.length; i++) {
            for (int j = i + 1; j < peaks.length; j++) {
                double mz1 = peaks[i].getMz();
                double mz2 = peaks[j].getMz();

                Double m1 = Adduct.getMonoisotopicMassFromMZ(mz1, "[M+H]+");
                Double m2 = Adduct.getMonoisotopicMassFromMZ(mz2, "[M+2H]2+");
                if (m1 != null && m2 != null && Math.abs(m1 - m2) <= toleranceDa) {
                    this.setAdduct("[M+H]+");
                    return;
                }

                Double m3 = Adduct.getMonoisotopicMassFromMZ(mz2, "[M+H]+");
                Double m4 = Adduct.getMonoisotopicMassFromMZ(mz1, "[M+2H]2+");
                if (m3 != null && m4 != null && Math.abs(m3 - m4) <= toleranceDa) {
                    this.setAdduct("[M+H]+");
                    return;
                }
            }
        }
    }

    private String findBestMatchingAdduct(double mz) {
        for (String adduct : AdductList.MAPMZPOSITIVEADDUCTS.keySet()) {
            Double mass = Adduct.getMonoisotopicMassFromMZ(mz, adduct);
            if (mass != null) return adduct;
        }
        return null;
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public IoniationMode getIonizationMode() {
        return ionizationMode;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    public double getNormalizedScore() {
        return (double) this.score / this.totalScoresApplied;
    }

    public int getCarbonCount() {
        return lipid.getCarbonCount();
    }

    public int getDoubleBondsCount() {
        return lipid.getDoubleBondsCount();
    }

    public String getLipidType() {
        return lipid.getLipidType();
    }

    public int getLipidClassRank() {
        switch (lipid.getLipidType()) {
            case "PG": return 1;
            case "PE": return 2;
            case "PI": return 3;
            case "PA": return 4;
            case "PS": return 5;
            case "PC": return 6;
            default: return 100;
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }
}