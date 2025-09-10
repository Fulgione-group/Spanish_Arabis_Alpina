/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Evolutionary_genetics;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author btjeng
 */
public class CalcDiv_DownSample_Allsites {

    public static void main(String[] args) throws IOException {
        System.out.println("Args: " + Arrays.toString(args));
//        System.out.println("0%");
//        calculate_diversity(new String[]{"--vcf", "/home/btjeng/Data/testing/test.vcf", "--cutoff", "0.1", "--popfile", "/home/btjeng/Data/testing/1000Genomes_final_Iberians.popfile", "--haploidize"});

        calculate_diversity(args);
    }
    //  }

    private static void calculate_diversity(String[] inputList) throws FileNotFoundException, IOException {
        String vcfpath = null;
        String outpath = null;
        String popfile = null;
        String replicate = "NA";
        String id = "NA";
        int seqlen = 0;
        int seqlen_bin = 0;
        Integer binsize = 10000;
        float cutoff = 0;
        boolean haploidize = false;
        Random r = new Random(System.currentTimeMillis());
        for (int i = 0; i < inputList.length; i++) {
//            System.out.println(inputList[i]);
            if (inputList[i].equals("--vcf") | inputList[i].equals("-i")) {
                vcfpath = inputList[(i + 1)];
//                System.out.println(vcfpath);
                i++;
            } else if (inputList[i].equals("--out") | inputList[i].equals("-o")) {
                outpath = inputList[(i + 1)];
                i++;
            } else if (inputList[i].equals("--cutoff") | inputList[i].equals("-c")) {
                cutoff = Float.parseFloat(inputList[(i + 1)]);
                i++;
            } else if (inputList[i].equals("--binsize") | inputList[i].equals("-s")) {
                binsize = Integer.parseInt(inputList[(i + 1)]);
                i++;
            } else if (inputList[i].equals("--haploidize") | inputList[i].equals("-h")) {
                haploidize = true;
            } else if (inputList[i].equals("--popfile") | inputList[i].equals("-p")) {
                popfile = inputList[(i + 1)];
                i++;
            }

        }
        if (outpath == null) {
            if (vcfpath.endsWith(".gz")) {
                outpath = vcfpath.substring(0, vcfpath.lastIndexOf(".")) + ".theta";
            } else {
                outpath = vcfpath + ".theta";
            }
        }
        if (vcfpath == null) {
            throw new IllegalArgumentException("No input specified");
        }
        BufferedWriter writer = null;
        BufferedWriter writerbin = null;
        writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath + ".gw.txt"))));
        writerbin = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath + ".bin.txt"))));
        InputStream inputStream = new BufferedInputStream(new FileInputStream(vcfpath));
        if (vcfpath.endsWith(".gz")) {
            inputStream = new GZIPInputStream(inputStream);
        }
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        BufferedReader reader2 = new BufferedReader(new InputStreamReader(new BufferedInputStream(new FileInputStream(popfile))));
        String[] linearray = null;
        String line;
        String chromosome = "start";
        int maxmissing = 0;
        int cap = 0;
        int ac = 0;
        int max_ac = 0;
        int bpos = 0;
        int endpos = 0;
        int startpos = 0;
        int[] maxmissing_array = null;
        int[] cap_array = null;
        int[] max_ac_array = null;
        double[] seqlen_array = null;
        int[] seqlen_bin_array = null;
        double[] pi_array = null;
        double[] pi_bin_array = null;
        double[] watterson_array = null;
        double[] watterson_bin_array = null;
        List<String> sampleIDs = new ArrayList<>();
        List<String> poptmplist = new ArrayList<>();
        List<Integer> popposlist = new ArrayList<>();
        List<List<Integer>> listOfpopposlist = new ArrayList<List<Integer>>();
        Map<String, List<String>> mapOfpoplist = new HashMap<>();
        Map<Integer, double[]> mapOfGWBin = new HashMap<>();
        while ((line = reader2.readLine()) != null) {
            linearray = P_TAB.split(line);
            if (linearray[1].equals("population")) {
                continue;
            } else {
                // if the population is already in the map add the sample to the list inside the map entry.
                if (mapOfpoplist.containsKey(linearray[1])) {
                    mapOfpoplist.get(linearray[1]).add(linearray[0]);
                } else {
                    //add new population to map and add to list inside the map entry.
                    poptmplist.add(linearray[0]);
                    mapOfpoplist.put(linearray[1], poptmplist);
                    poptmplist = new ArrayList<>();
                }
            }
        }
        writerbin.write("chromosome" + "\t" + "start" + "\t" + "end" + "\t" + "population" + "\t" + "number_of_genotypes" + "\t" + "number_of_nonmissingsites" + "\t" + "segragating_sites" + "\t" + "watterson" + "\t" + "pi" + "\t" + "TajimasD" + "\t");
        writerbin.newLine();
        while ((line = reader.readLine()) != null) {
            if (line.charAt(0) != '#') {
                linearray = P_TAB.split(line);
                int newstart = 0;
                while (!chromosome.equals(linearray[0]) | Integer.parseInt(linearray[1]) > (startpos + binsize)) {
                    //just to bypas this code for the first position.
                    if (!chromosome.equals("start")) {
                        if (chromosome.equals(linearray[0])) {
                            endpos = startpos + binsize - 1;
                            newstart = startpos + binsize;
                        } else {
                            endpos = bpos;
                            newstart = 1;
                        }

                        // loop over all populations
                        for (int k = 0; k < listOfpopposlist.size(); k++) {
                            double nsegsites = watterson_bin_array[k];
                            double watterson = watterson_bin_array[k] / harmonicNumber(max_ac_array[k]-1);
                            double tajV = tajDvar(max_ac_array[k], nsegsites);
                            double tajD = (pi_bin_array[k] - watterson) / Math.sqrt(tajV);
//                            if(k==(listOfpopposlist.size()-1)){
//                                System.out.println(max_ac_array[k] + "\t" + nsegsites);
//                                System.out.println(tajDvar(5, 4));
//                            }
                            watterson = watterson / seqlen_bin_array[k];
                            double pi = pi_bin_array[k] / seqlen_bin_array[k];
                            writerbin.write(chromosome + "\t" + startpos + "\t" + endpos + "\t" + mapOfpoplist.keySet().toArray()[k] + "\t" + max_ac_array[k] + "\t" +  seqlen_bin_array[k] + "\t" + nsegsites + "\t" + watterson + "\t" + pi + "\t" + tajD);
                            writerbin.newLine();
                        }
                        //creat new map of frequencies for a bin
                        Arrays.fill(pi_bin_array, 0);
                        Arrays.fill(watterson_bin_array, 0);
                        Arrays.fill(seqlen_bin_array, 0);
                    } else {
                        //this will let the other lines pass
                        newstart = 1;
                    }                        //new start position is the new 
                    startpos = newstart;
                    chromosome = linearray[0];
                    bpos = Integer.parseInt(linearray[1]);
                }
                if (Integer.parseInt(linearray[1]) <= (startpos + binsize)) {
                    ac = 0;
                    boolean missingpass = false;
                    int[] missing_array = new int[listOfpopposlist.size()];
                    int[] ac_array = new int[listOfpopposlist.size()];
                    Arrays.fill(missing_array, 0);
                    List<List<String>> Nonmissing_list = new ArrayList<>();
                    for (int k = 0; k < listOfpopposlist.size(); k++) {
                        Nonmissing_list.add(new ArrayList<String>());
                    }
                    for (int k = 0; k < listOfpopposlist.size(); k++) {
                        for (int i = 0; i < listOfpopposlist.get(k).size(); i++) {
                            if (linearray[listOfpopposlist.get(k).get(i)].charAt(0) == '.') {
                                missing_array[k]++;
                            } else {
//                                System.out.println(linearray[listOfpopposlist.get(k).get(i)]);
                                Nonmissing_list.get(k).add(linearray[listOfpopposlist.get(k).get(i)]);
                            }
                        }
                    }
                    for (int k = 0; k < listOfpopposlist.size(); k++) {
                        if (missing_array[k] > maxmissing_array[k]) {
                            continue;
                        } else {
//                            System.out.println("yeah");
                            seqlen_array[k]++;
                            seqlen_bin_array[k]++;
//                            System.out.println(Nonmissing_list.get(k).size());
                            reduceArray(Nonmissing_list.get(k), r, cap_array[k]);
//                    System.out.println(Nonmissing_list.size());
                            if (haploidize) {
//                                System.out.println(Nonmissing_list.get(k).size());
//                                System.out.println(Nonmissing_list.size());
                                for (int i = 0; i < Nonmissing_list.get(k).size(); i++) {
                                    int randomOfTwoInts = r.nextBoolean() ? 0 : 2;
//                                    System.out.println(Nonmissing_list.get(k).get(i).charAt(randomOfTwoInts));
                                    if (Nonmissing_list.get(k).get(i).charAt(randomOfTwoInts) == '1') {
//                                        System.out.println(Nonmissing_list.get(k).get(i).charAt(randomOfTwoInts));
                                        ac_array[k]++;
                                    }
                                }
//                                System.out.println(ac_array[k]);
                            } else {
                                for (int i = 0; i < Nonmissing_list.get(k).size(); i++) {
//                        System.out.println(Nonmissing_list.get(i)+ "\t" + Nonmissing_list.get(i).charAt(0) + "\t" +Nonmissing_list.get(i).charAt(2));
                                    if (Nonmissing_list.get(k).get(i).charAt(0) == '1') {
                                        ac_array[k]++;
//                            System.out.println("first alt");
                                    }
                                    if (Nonmissing_list.get(k).get(i).charAt(2) == '1') {
                                        ac_array[k]++;
//                            System.out.println("second alt");
                                    }
                                }
                            }
                            if (ac_array[k] > 0 & ac_array[k] < max_ac_array[k]) {
//                                System.out.println("yes SNP");
                                pi_array[k] += (double) 2 * ac_array[k] * (max_ac_array[k] - ac_array[k]) / (max_ac_array[k] * (max_ac_array[k] - 1));
                                pi_bin_array[k] += (double) 2 * ac_array[k] * (max_ac_array[k] - ac_array[k]) / (max_ac_array[k] * (max_ac_array[k] - 1));
                                watterson_array[k]++;
                                watterson_bin_array[k]++;
                            }
                        }
                    }
                }
                bpos = Integer.parseInt(linearray[1]);
            } else if (line.startsWith("#CHROM")) {
                linearray = P_TAB.split(line);
                for (String k : mapOfpoplist.keySet()) {
                    for (int j = 0; j < mapOfpoplist.get(k).size(); j++) {
                        for (int i = 0; i < linearray.length; i++) {
                            if (mapOfpoplist.get(k).get(j).equals(linearray[i])) {
                                //checked this part
//                                    if(k.equals("NO06")){
//                                        System.out.println("colnum: " +i);
//                                    }
                                popposlist.add(i);
                            }
                        }
                    }
                    listOfpopposlist.add(popposlist);
                    popposlist = new ArrayList<>();
                }
                maxmissing_array = new int[listOfpopposlist.size()];
                cap_array = new int[listOfpopposlist.size()];
                max_ac_array = new int[listOfpopposlist.size()];
                seqlen_array = new double[listOfpopposlist.size()];
                seqlen_bin_array = new int[listOfpopposlist.size()];
                pi_array = new double[listOfpopposlist.size()];
                watterson_array = new double[listOfpopposlist.size()];
                pi_bin_array = new double[listOfpopposlist.size()];
                watterson_bin_array = new double[listOfpopposlist.size()];
                Arrays.fill(seqlen_array, 0);
                Arrays.fill(seqlen_bin_array, 0);
                Arrays.fill(pi_array, 0);
                Arrays.fill(watterson_array, 0);
                Arrays.fill(pi_bin_array, 0);
                Arrays.fill(watterson_bin_array, 0);
                for (int k = 0; k < listOfpopposlist.size(); k++) {
                    maxmissing_array[k] = (int) (cutoff * (float) listOfpopposlist.get(k).size() + 0.5);
                    cap_array[k] = listOfpopposlist.get(k).size() - maxmissing_array[k];
//                        System.out.println(listOfpopposlist.get(k).size() + "\t"+ maxmissing +"\t"  +maxmissingarray[k]);
                    if (haploidize) {
                        max_ac_array[k] = cap_array[k];
                    } else {
                        max_ac_array[k] = 2 * cap_array[k];
                    }
                }
//                System.out.println("cap" + "\t" + cap + "\t" + "maxmissing" +"\t" + maxmissing + "\t" +"linelength" + "\t" + linearray.length);
            }
        }
        endpos = bpos;
        for (int k = 0; k < listOfpopposlist.size(); k++) {
            double nsegsites = watterson_bin_array[k];
            double watterson = watterson_bin_array[k] / harmonicNumber(max_ac_array[k]-1);
            double tajV = tajDvar(max_ac_array[k], nsegsites);
            double tajD = (pi_bin_array[k] - watterson) / Math.sqrt(tajV);
            watterson = watterson / seqlen_bin_array[k];
            double pi = pi_bin_array[k] / seqlen_bin_array[k];
            writerbin.write(chromosome + "\t" + startpos + "\t" + endpos + "\t" + mapOfpoplist.keySet().toArray()[k] + "\t" + max_ac_array[k] + "\t" +  seqlen_bin_array[k]+ "\t" + nsegsites + "\t" + watterson + "\t" + pi + "\t" + tajD);
            writerbin.newLine();
        }
        writer.write("population" + "\t" + "number_of_genotypes" + "\t" + "number_of_nonmissingsites" +"\t" + "segragating_sites" + "\t" + "watterson" + "\t" + "pi" + "\t" + "TajimasD" + "\t");
        writer.newLine();
        for (int k = 0; k < listOfpopposlist.size(); k++) {
            double nsegsites = watterson_array[k];
            double watterson = watterson_array[k] / harmonicNumber(max_ac_array[k]-1);
            double tajV = tajDvar(max_ac_array[k], nsegsites);
            double tajD = (pi_array[k] - watterson) / Math.sqrt(tajV);
            watterson = watterson / seqlen_array[k];
            double pi = pi_array[k] / seqlen_array[k];
            writer.write(mapOfpoplist.keySet().toArray()[k] + "\t" + max_ac_array[k] + "\t" +  seqlen_array[k] + "\t" + nsegsites + "\t" + watterson + "\t" + pi + "\t" + tajD);
            writer.newLine();
        }
        writerbin.flush();
        writerbin.close();
        writer.flush();
        writer.close();
    }

    private static double harmonicNumber(int n) {
        double d = n;
        // Base Cases
        if (d == 0) {
            return 0;
        }
        // Recur
        return (1 / d) + harmonicNumber(n - 1);
    }

    private static double harmonicNumberSq(int n) {
        double d = n * n;
        // Base Cases
        if (d == 0) {
            return 0;
        }
        // Recur
        return (1 / d) + harmonicNumberSq(n - 1);
    }

    private static void reduceArray(List<String> list, Random r, int cap) {
        while (list.size() > cap) {
            int index = r.nextInt(list.size());
            list.remove(index);
        }
    }

    private static double tajDvar(int ns, double s) {
        double n = ns;
        double a1 = harmonicNumber((ns - 1));
        double a2 = harmonicNumberSq((ns - 1));
        double b1 = (n + 1) / (3 * (n - 1));
        double b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1));
        double c1 = b1 - 1 / a1;
        double c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1);
        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);
        double var = e1 * s + e2 * s * (s - 1);
        return (var);
    }

    private static final Pattern P_TAB = Pattern.compile("\t");
    private static final Pattern P_SPACE = Pattern.compile(" ");
}
