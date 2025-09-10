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
public class CalcDiv_and_Fis_DownSample_Allsites {

    public static void main(String[] args) throws IOException {
        System.out.println("Args: " + Arrays.toString(args));
//        System.out.println("0%");
//        calculate_diversity(new String[]{"--vcf", "/home/btjeng/Data/testing/test2.vcf", "--cutoff", "0.1", "--popfile", "/home/btjeng/Data/testing/unrelated_pop_fiile.txt"});
//        calculate_diversity(new String[]{"--vcf", "/home/btjeng/Data/testing/test2.vcf", "--cutoff", "0.1", "--popfile", "/home/btjeng/Data/testing/unrelated_pop_file_ES0304.txt"});
        calculate_diversity(args);
    }
    //  }

    private static void calculate_diversity(String[] inputList) throws FileNotFoundException, IOException {
        String vcfpath = null;
        String outpath = null;
        String popfile = null;
        Integer binsize = 10000;
        float cutoff = 0;
        Random r = new Random(System.currentTimeMillis());
        for (int i = 0; i < inputList.length; i++) {
//            System.out.println(inputList[i]);
            if (inputList[i].equals("--vcf") | inputList[i].equals("-i")) {
                vcfpath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--out") | inputList[i].equals("-o")) {
                outpath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--cutoff") | inputList[i].equals("-c")) {
                cutoff = Float.parseFloat(inputList[(i + 1)]);
                i++;
            }
            if (inputList[i].equals("--binsize") | inputList[i].equals("-s")) {
                binsize = Integer.parseInt(inputList[(i + 1)]);
                i++;
            }
            if (inputList[i].equals("--popfile") | inputList[i].equals("-p")) {
                popfile = inputList[(i + 1)];
                i++;
            }
        }
        if (outpath == null) {
            if (vcfpath.endsWith(".gz")) {
                outpath = vcfpath.substring(0, vcfpath.lastIndexOf(".")) + ".SummaryStats";
            } else {
                outpath = vcfpath + ".SummaryStats";
            }
        }
        if (vcfpath == null) {
            throw new IllegalArgumentException("No input specified");
        }
        BufferedWriter writer = null;
        BufferedWriter writer2 = null;
        BufferedWriter writer3 = null;
        writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath + ".diversity.txt"))));
        writer2 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath + ".Fstat.txt"))));
        writer3 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath + ".bin.Fstat.txt"))));
        InputStream inputStream = new BufferedInputStream(new FileInputStream(vcfpath));
        if (vcfpath.endsWith(".gz")) {
            inputStream = new GZIPInputStream(inputStream);
        }
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        BufferedReader reader2 = new BufferedReader(new InputStreamReader(new BufferedInputStream(new FileInputStream(popfile))));
        String[] linearray = null;
        String line;
        Map<String, List<String>> mapOfpoplist = new HashMap<>();
        List<String> poptmplist = new ArrayList<>();
        List<Integer> popposlist = new ArrayList<>();
        List<List<Integer>> listOfpopposlist = new ArrayList<List<Integer>>();
        Map<Integer, double[]> mapOfGWBin = new HashMap<>();
        List<String> sampleIDs = new ArrayList<>();
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
        List<Integer> maxmissing = new ArrayList<>();
        List<Integer> cap = new ArrayList<>();
        List<String> pops_present = new ArrayList<>();
        Double[] nsegsites = null;
        Double[] watterson = null;
        Double[] pi = null;
        Double[] tajV = null;
        Double[] tajD = null;
        Double[] seqlen = null;
        Double ho_tmp = (double) 0;
        Double ho = (double) 0;
        Double hs = (double) 0;
        Double hs_tmp = (double) 0;
        Double ht = (double) 0;
        Double ht_tmp = (double) 0;
        Double dst = (double) 0;
        int ac = 0;
        Double hom = null;
        Double af = null;
        Double hom_sum = (double) 0;
        Double pop_sum = (double) 0;
        Double sum_prop = (double) 0;
        Double n_nei = (double) 0;
        Double p_nei = (double) 0;
        Double p_nei_r = (double) 0;
        Double p_nei_a = (double) 0;
        Double p2_nei = (double) 0;
        Double p2_nei_r = (double) 0;
        Double p2_nei_a = (double) 0;
        Double Fis = (double) 0;
        Double Fst = (double) 0;
        Double SNP_counter = (double) 0;
        Double ho_bin = (double) 0;
        Double hs_bin = (double) 0;
        Double ht_bin = (double) 0;
        Double dst_bin = (double) 0;
        Double Fis_bin = (double) 0;
        Double Fst_bin = (double) 0;
        Integer bincounter = 0;
        writer.write("populations" + "\t" + "number_of_genotypes" + "\t" + "segragating_sites" + "\t" + "watterson" + "\t" + "pi" + "\t" + "TajimasD");
        writer.newLine();
        writer3.write("bincount" + "\t" + "Ho" + "\t" + "Hs" + "\t" + "Ht" + "\t" + "Dst" + "\t" + "Fis" + "\t" + "Fst");
        writer3.newLine();
        writer2.write("Ho" + "\t" + "Hs" + "\t" + "Ht" + "\t" + "Dst" + "\t" + "Fis" + "\t" + "Fst");
        writer2.newLine();
        while ((line = reader.readLine()) != null) {
            linearray = P_TAB.split(line);
            if (line.startsWith("#CHROM")) {
                boolean contains_pop = false;
                for (String k : mapOfpoplist.keySet()) {
                    contains_pop = false;
                    for (int j = 0; j < mapOfpoplist.get(k).size(); j++) {
                        for (int i = 9; i < linearray.length; i++) {
                            if (mapOfpoplist.get(k).get(j).equals(linearray[i])) {
                                popposlist.add(i);
                                contains_pop = true;
                            }
                        }
                    }
                    if (contains_pop) {
//                        System.out.println(k);
                        pops_present.add(k);
                    }
                    listOfpopposlist.add(popposlist);
                    maxmissing.add((int) (cutoff * popposlist.size() + 0.5));
                    cap.add(popposlist.size() - (int) (cutoff * popposlist.size() + 0.5));
                    popposlist = new ArrayList<>();
                }
                nsegsites = new Double[listOfpopposlist.size()];
                watterson = new Double[listOfpopposlist.size()];
                pi = new Double[listOfpopposlist.size()];
                tajV = new Double[listOfpopposlist.size()];
                tajD = new Double[listOfpopposlist.size()];
                seqlen = new Double[listOfpopposlist.size()];
                Arrays.fill(nsegsites, (double) 0);
                Arrays.fill(watterson, (double) 0);
                Arrays.fill(pi, (double) 0);
                Arrays.fill(seqlen, (double) 0);
                Arrays.fill(tajV, (double) 0);
                Arrays.fill(tajD, (double) 0);
            } else if (line.charAt(0) != '#') {
                if (linearray[4].contains(",")) {
//                    System.out.println(linearray[4]);
                    continue;
                }
                p_nei = (double) 0;
                p_nei_r = (double) 0;
                p_nei_a = (double) 0;
                p2_nei = (double) 0;
                p2_nei_r = (double) 0;
                p2_nei_a = (double) 0;
                n_nei = (double) 0;
                ho_tmp = (double) 0;
                hs_tmp = (double) 0;
                ht_tmp = (double) 0;
                hom_sum = (double) 0;
                pop_sum = (double) 0;
                sum_prop = (double) 0;
                boolean countsite = false;
                for (int i = 0; i < listOfpopposlist.size(); i++) {
                    String[] subarray = new String[listOfpopposlist.get(i).size()];
                    for (int j = 0; j < listOfpopposlist.get(i).size(); j++) {
                        subarray[j] = linearray[listOfpopposlist.get(i).get(j)];
                    }
//                    System.out.println(Arrays.toString(subarray) + "\t" + i);
                    Double[] CurrentStats = Stats(subarray, maxmissing.get(i), cap.get(i));
                    seqlen[i] += CurrentStats[0];
                    pi[i] += CurrentStats[1];
                    nsegsites[i] += CurrentStats[2];
                    hom = CurrentStats[3];
                    af = CurrentStats[4];
                    if (hom >= 0) {
                        hom_sum += hom;
                        pop_sum++;
                        sum_prop += (double) 1 / cap.get(i);
//                        System.out.println(sum_prop + "\t" + cap.get(i));
                        p_nei_a += af;
                        p_nei_r += (1 - af);
                        p2_nei_a += af * af;
                        p2_nei_r += (1 - af) * (1 - af);
                        countsite = true;
                        SNP_counter++;
//                        System.out.println(hom + "\t"+ af);
                    }
                }
                // weighted stats from nei for Fis and Fst per loci later ratio of sums.
                if (countsite) {
                    p_nei_a = p_nei_a / pop_sum;
                    p_nei_r = p_nei_r / pop_sum;
                    p2_nei = (p2_nei_r + p2_nei_a) / pop_sum;
//                System.out.println(p2_nei);
                    n_nei = pop_sum / sum_prop;
//                System.out.println(n_nei + "\t" + sum_prop);
                    ho_tmp = 1 - hom_sum / pop_sum;
                    ho += ho_tmp;
                    ho_bin += ho_tmp;
                    hs_tmp = (n_nei / (n_nei - 1)) * (1 - p2_nei - ho_tmp / (2 * n_nei));
                    hs += hs_tmp;
                    hs_bin += hs_tmp;
                    ht_tmp = 1 - (p_nei_a * p_nei_a + p_nei_r * p_nei_r) + (hs_tmp / (n_nei * pop_sum)) - (ho_tmp / (2 * n_nei * pop_sum));
                    ht += ht_tmp;
                    ht_bin += ht_tmp;
                    dst += ht_tmp - hs_tmp;
                    dst_bin += ht_tmp - hs_tmp;
//                System.out.println(1-p2_nei +"\t" + ho_tmp);
//                    System.out.println(ht_tmp + "\t" + hs_tmp + "\t");
                }
                if (SNP_counter % binsize == 0) {
                    bincounter++;
                    Fis_bin = 1 - ho_bin / hs_bin;
                    Fst_bin = dst_bin / ht_bin;
                    writer3.write(bincounter + "\t" + ho_bin + "\t" + hs_bin + "\t" + ht_bin + "\t" + dst_bin + "\t" + Fis_bin + "\t" + Fst_bin + "\t");
                    writer3.newLine();
                    ho_bin = (double) 0;
                    hs_bin = (double) 0;
                    ht_bin = (double) 0;
                    dst_bin = (double) 0;
                }
            }
        }
//        dst = ht - hs;
//        System.out.println(dst + "\t" + ht);
        Fis = 1 - ho / hs;
        Fst = dst / ht;
        writer2.write(ho + "\t" + hs + "\t" + ht + "\t" + dst + "\t" + Fis + "\t" + Fst + "\t");
        writer2.newLine();
        writer2.flush();
        writer2.close();
        for (int i = 0; i < listOfpopposlist.size(); i++) {
            watterson[i] = nsegsites[i] / harmonicNumber(2 * cap.get(i));
            tajV[i] = tajDvar(2 * cap.get(i), nsegsites[i]);
            tajD[i] = (pi[i] - watterson[i]) / Math.sqrt(tajV[i]);
            watterson[i] = watterson[i] / seqlen[i];
            pi[i] = pi[i] / seqlen[i];
            writer.write(pops_present.get(i) + "\t" + cap.get(i) + "\t" + nsegsites[i] + "\t" + watterson[i] + "\t" + pi[i] + "\t" + tajD[i]);
            writer.newLine();
        }
        writer.flush();
        writer.close();
        bincounter++;
        Fis_bin = 1 - ho_bin / hs_bin;
        Fst_bin = dst_bin / ht_bin;
        writer3.write(bincounter + "\t" + ho_bin + "\t" + hs_bin + "\t" + ht_bin + "\t" + dst_bin + "\t" + Fis_bin + "\t" + Fst_bin + "\t");
        writer3.newLine();
        ho_bin = (double) 0;
        hs_bin = (double) 0;
        ht_bin = (double) 0;
        dst_bin = (double) 0;
        writer3.flush();
        writer3.close();
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

    private static Double[] Stats(String[] subarray, Integer maxmissing_tmp, Integer cap_tmp) {
        Random r = new Random(System.currentTimeMillis());
        double ac = 0;
        boolean missingpass = false;
        double pi_tmp = 0;
        double watterson_tmp = 0;
        double missing = 0;
        double hom = -1;
        double seqlen_tmp = 0;
        List<String> Nonmissing_list = new ArrayList<>();
        for (int i = 0; i < subarray.length; i++) {
            if (subarray[i].charAt(0) == '.') {
                missing++;
            } else {
                Nonmissing_list.add(subarray[i]);
            }
        }
        if (missing > maxmissing_tmp) {
            missingpass = false;
//                    System.out.println(seqlen + "\t" +missing);
//                    System.out.println(missing + "\t" +line);
        } else {
            seqlen_tmp++;
            reduceArray(Nonmissing_list, r, cap_tmp);
            hom = 0;
//                    System.out.println(Nonmissing_list.size());
            for (int i = 0; i < Nonmissing_list.size(); i++) {
//                        System.out.println(Nonmissing_list.get(i)+ "\t" + Nonmissing_list.get(i).charAt(0) + "\t" +Nonmissing_list.get(i).charAt(2));
                if (Nonmissing_list.get(i).charAt(0) == '1') {
                    ac++;
//                            System.out.println("first alt");
                }
                if (Nonmissing_list.get(i).charAt(2) == '1') {
                    ac++;
//                            System.out.println("second alt");
                }
                if (Nonmissing_list.get(i).charAt(2) == Nonmissing_list.get(i).charAt(0)) {
                    hom++;
//                    System.out.println(Nonmissing_list.get(i).charAt(2) + "\t" +Nonmissing_list.get(i).charAt(0));
                }
//                else{
//                      System.out.println(Nonmissing_list.get(i).charAt(2) + "\t" +Nonmissing_list.get(i).charAt(0));
//                }
            }
//                    System.out.println(2*cap);
            if (ac > 0 & ac < (2 * cap_tmp)) {
                pi_tmp = (double) 2 * ac * (2 * cap_tmp - ac) / (2 * cap_tmp * (2 * cap_tmp - 1));
//                        System.out.println( cap + "\t" + ac + "\t" +(double) 2 * ac * (2 * cap - ac) / (2* cap * (2* cap- 1)));
//                        System.out.println(ac);
                watterson_tmp++;
            }
//            if(hom==0){
//                System.out.println(Nonmissing_list.toString());
//            }
//            System.out.println(hom + "\t" + cap_tmp);
            hom = hom / cap_tmp;
//            System.out.println(hom);
        }
        return new Double[]{seqlen_tmp, pi_tmp, watterson_tmp, hom, ac / (2 * cap_tmp)};
    }

    private static final Pattern P_TAB = Pattern.compile("\t");
    private static final Pattern P_SPACE = Pattern.compile(" ");
}
