package adduct;

public class Adduct {
    private static int extractMultimer(String adduct){
        java.util.regex.Matcher matcher = java.util.regex.Pattern.compile("\\[(\\d*)M").matcher(adduct);
        if (matcher.find()){
            String multimerStr= matcher.group(1);
            return multimerStr.isEmpty() ?1: Integer.parseInt(multimerStr);
        }
        return 1;
    }
    private static int extractCharge(String adduct){
        java.util.regex.Matcher matcher = java.util.regex.Pattern.compile("(\\d*)([+-])\\]").matcher(adduct);
        if (matcher.find()) {
            String chargeStr = matcher.group(1);
            return (chargeStr == null || chargeStr.isEmpty()) ? 1 : Integer.parseInt(chargeStr);
        }
        return 1;
    }

    /**
     * Calculate the mass to search depending on the adduct hypothesis
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, String adduct) {

        Double adductMass = AdductList.MAPMZPOSITIVEADDUCTS.getOrDefault(adduct, 0.0);
        if (adductMass == 0.0) return null;


        int multimer = extractMultimer(adduct);


        int charge = extractCharge(adduct);


        if (charge == 0 || multimer == 0) return null;


        return (mz + adductMass) * charge / multimer;
    }


    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     *
     * @param monoisotopicMass
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, String adduct) {

        Double adductMass = AdductList.MAPMZPOSITIVEADDUCTS.getOrDefault(adduct, 0.0);
        if (adductMass == 0.0) return null;


        int multimer = extractMultimer(adduct);
        int charge = extractCharge(adduct);

        // Validación
        if (charge == 0 || multimer == 0) return null;


        return (monoisotopicMass * multimer / charge) - adductMass;
    }


    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param theoreticalMass Theoretical mass of the compound
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass) * 1000000
                / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param measuredMass    Mass measured by MS
     * @param ppm ppm of tolerance
     */
    public static double calculateDeltaPPM(Double experimentalMass, int ppm) {
        double deltaPPM;
        deltaPPM =  Math.round(Math.abs((experimentalMass * ppm) / 1000000));
        return deltaPPM;

    }




}
