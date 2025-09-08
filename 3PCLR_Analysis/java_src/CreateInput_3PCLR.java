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
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author btjeng
 */
public class CreateInput_3PCLR {

    public static void main(String[] args) throws FileNotFoundException {
        System.out.println("Args: " + Arrays.toString(args));
//        downsampledAC(new String[]{"--haploidize", "--map", "/home/btjeng/Data/Spanish_adaptation/geneticMap_chr1_3pclr_rates.txt", "--vcf", "/home/btjeng/Data/Spanish_adaptation/Set1_test.vcf", "--popfile", "/home/btjeng/Data/Spanish_adaptation/1000Genomes_final_Canabria_3pxpclr.popfile", "--listofpop", "ES03,ES04,FRANCE", "--cutoff", "0.1", "--chr", "chr1"});
        downsampledAC(args);
    }

    // vcfpath1 = full path to the vcf
    // popfile = full path to popfile specifying a population ID for every sample ID separated by spaces. Each sample is writen on one line.
    // pop1 = population ID of population1
    // pop2 = population ID of population2
    // ref = sample ID of the outgroup
    // output = full path of the outputfile including its name
    private static void downsampledAC(String[] inputList) throws FileNotFoundException {
        String inpath = null;
        String mappath = null;
        String popfile = null;
        String listofpop = null;
        String chromosome = null;
        String listofpopout = null;
        String ref = null;
        String outpath = null;
        String[] arrayofpop = null;
        boolean haploidize = false;
        float cutoff = 0;
        List<String> sampleIDs = new ArrayList<>();
        int refpos = 0;
        Random r = new Random(System.currentTimeMillis());
        for (int i = 0; i < inputList.length; i++) {
            // full input path of the vcf file
            if (inputList[i].equals("--vcf") | inputList[i].equals("-i")) {
                inpath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--map") | inputList[i].equals("-m")) {
                mappath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--popfile") | inputList[i].equals("-pf")) {
                popfile = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--listofpop") | inputList[i].equals("-lp")) {
                listofpop = inputList[(i + 1)];
                listofpopout = listofpop.replace(",", ".");
                i++;
            }
            if (inputList[i].equals("--ref") | inputList[i].equals("-r")) {
                ref = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--cutoff") | inputList[i].equals("-c")) {
                cutoff = Float.parseFloat(inputList[(i + 1)]);
                i++;
            }
            //full path of the output vcf 
            if (inputList[i].equals("--out") | inputList[i].equals("-o")) {
                outpath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--chr") | inputList[i].equals("-c")) {
                chromosome = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--haploidize") | inputList[i].equals("-h")) {
                haploidize = true;
            }
        }
        if (outpath == null) {
            if (inpath.endsWith(".gz")) {
                outpath = inpath.substring(0, inpath.lastIndexOf(".")) + "." + listofpopout +  "." + chromosome + "."+ ".3PCLRinput";
            } else {
                outpath = inpath + "." + listofpopout  + "." + chromosome + "."+ ".3PCLRinput";
            }
        }
        if (inpath == null) {
            throw new IllegalArgumentException("No input specified");
        }
        arrayofpop = P_COMMA.split(listofpop);
        ArrayList<ArrayList<String>> popsamplelist = new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<Integer>> popsampleposlist = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> popsamplepostmp = new ArrayList<ArrayList<Integer>>();
        Integer[] SFS = null;
        Integer[] maxmissing = null;
        Integer[] missing = null;
        Integer[] acount = null;
        Integer[] dim = null;
        String[] linearray = null;
        String[] maparray = null;
        String refsplit = null;
        String line;
        String map_line;
        double rrate = 0;
        double ppos = 0;
        double gpos = 0;
        double gpossite = 0;
        boolean isMultiallelicOrInvariable = true;
        for (int i = 0; i < arrayofpop.length; i++) {
            popsamplelist.add(new ArrayList<String>());
            popsampleposlist.add(new ArrayList<Integer>());
            popsamplepostmp.add(new ArrayList<Integer>());
        }
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath))));
        try ( BufferedReader reader2 = new BufferedReader(new FileReader(popfile))) {
            while ((line = reader2.readLine()) != null) {
                linearray = P_SPACE.split(line);
//                    linearray = line.split(" ");
//                    System.out.println(Arrays.toString(linearray));
//                    System.out.println(linearray[0]);
//                    System.out.println(linearray[1]);
                for (int i = 0; i < arrayofpop.length; i++) {
                    if (linearray[1].equals(arrayofpop[i])) {
                        popsamplelist.get(i).add(linearray[0]);
                    }
                }
            }
//            for(int i = 0 ; i< popsamplelist.size(); i++){
//                System.out.println(popsamplelist.get(i).size());
//            }
            maxmissing = new Integer[popsamplelist.size()];
            missing = new Integer[popsamplelist.size()];
            acount = new Integer[popsamplelist.size()];
        } catch (IOException ex) {
            Logger.getLogger(Main.class
                    .getName()).log(Level.SEVERE, null, ex);
        }

        String genoarray = null;
        //go through every line of the vcf
        try {
            InputStream inputStream = new BufferedInputStream(new FileInputStream(inpath));
            if (inpath.endsWith(".gz")) {
                inputStream = new GZIPInputStream(inputStream);
            }
            BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
            InputStream inputStream3 = new BufferedInputStream(new FileInputStream(mappath));
            BufferedReader reader3 = new BufferedReader(new InputStreamReader(inputStream3));
            int helparray[] = {0, 2};
            map_line = reader3.readLine();
            maparray = P_TAB.split(map_line);
            rrate = Double.parseDouble(maparray[3]);
            ppos = Double.parseDouble(maparray[1]);
            gpos = Double.parseDouble(maparray[2]);
            while ((line = reader.readLine()) != null) {
                List<String> pop1tmp = new ArrayList<>();
                List<String> pop2tmp = new ArrayList<>();
//                System.out.println(chromosome);
                //Disregard header lines
                if ((line.charAt(0) != '#')) {
                    isMultiallelicOrInvariable = linearray[4].contains(",") || linearray[4].equals(".");
                    linearray = P_TAB.split(line);
                    if (isMultiallelicOrInvariable || !linearray[0].equals(chromosome)) {
                        continue;
                    }
//                    System.out.println(Double.parseDouble(linearray[1]) + "\t" + ppos + "\t" + rrate);
                    while (Double.parseDouble(linearray[1]) > ppos && (map_line = reader3.readLine()) != null) {
//                        System.out.println(Double.parseDouble(linearray[1]) > ppos );
                        maparray = P_TAB.split(map_line);
                        rrate = Double.parseDouble(maparray[3]);
                        ppos = Double.parseDouble(maparray[1]);
                        gpos = Double.parseDouble(maparray[2]);
                    }
//                    System.out.println(Double.parseDouble(linearray[1]) + "\t" + ppos);
                    gpossite = gpos - rrate * (ppos - Double.parseDouble(linearray[1]));
                    if (ref != null) {
                        refsplit = linearray[refpos];
                        refsplit.charAt(0);
                        refsplit.charAt(2);
                    } else {
                        refsplit = "0/0";
                    }
                    //calculate missingness of a SNP
                    boolean missingpass = true;
                    Arrays.fill(missing, 0);
                    if (ref != null) {
                        if (refsplit.charAt(0) == '.' || refsplit.charAt(0) != refsplit.charAt(2)) {
                            continue;
                        }
                    }
                    for (int i = 0; i < popsampleposlist.size(); i++) {
                        popsamplepostmp.get(i).clear();  // Clear temporary list before filling it again
                        for (int j = 0; j < popsampleposlist.get(i).size(); j++) {
                            if (linearray[popsampleposlist.get(i).get(j)].charAt(0) == '.') {
                                missing[i]++;
                            } else {
                                popsamplepostmp.get(i).add(popsampleposlist.get(i).get(j));
                            }
                        }
                        if (missing[i] > maxmissing[i]) {
                            missingpass = false;
                        }
                    }
                    //if the SNP has enough genotypes called proceed and randomly downsample
                    if (missingpass) {
                        Arrays.fill(acount, 0);
//                        int acount = 0;
//                        int acount2 = 0;
                        for (int i = 0; i < popsamplepostmp.size(); i++) {
                            reduceArray(popsamplepostmp.get(i), r, (popsampleposlist.get(i).size() - maxmissing[i]));
                            for (int j = 0; j < popsamplepostmp.get(i).size(); j++) {
                                genoarray = linearray[popsamplepostmp.get(i).get(j)];
                                if (haploidize) {
                                    int randomOfTwoInts = r.nextBoolean() ? 0 : 2;
//                                    System.out.println(randomOfTwoInts);
                                    if (genoarray.charAt(randomOfTwoInts) != refsplit.charAt(0)) {
                                        acount[i]++;
                                    }
                                } else {
                                    if (genoarray.charAt(0) != refsplit.charAt(0)) {
                                        acount[i]++;
                                    }
                                    if (genoarray.charAt(2) != refsplit.charAt(0)) {
                                        acount[i]++;
                                    }
                                }
                            }
                        }
                        if (acount[2] == 0) {
                            continue;
                        }
                        writer.write(linearray[0] + "\t" + linearray[1] + "\t" + gpossite + "\t" + acount[0] + "\t" + (popsampleposlist.get(0).size() - maxmissing[0]) + "\t" + acount[1] + "\t" + (popsampleposlist.get(1).size() - maxmissing[1]) + "\t" + acount[2] + "\t" + (popsampleposlist.get(2).size() - maxmissing[2]));
                        writer.newLine();
                    }
                    //go through the vcf header line that contains the sample IDs
                } else if (line.startsWith("#CHROM")) {
                    linearray = P_TAB.split(line);
                    for (int i = 0; i < linearray.length; i++) {
                        for (int j = 0; j < popsamplelist.size(); j++) {
                            for (int k = 0; k < popsamplelist.get(j).size(); k++) {
                                if (popsamplelist.get(j).get(k).equals(linearray[i])) {
                                    popsampleposlist.get(j).add(i);
                                } else if (ref != null && linearray[i].equals(ref)) {
                                    refpos = i;
                                }
                            }
                        }
                    }
                    ArrayList<Integer> sizes = new ArrayList<Integer>();
                    BigInteger arraysize = BigInteger.valueOf(1);
                    for (int i = 0; i < popsampleposlist.size(); i++) {
                        maxmissing[i] = (int) Math.floor(popsampleposlist.get(i).size() * cutoff);
//                        System.out.println(popsampleposlist.get(i).size());
                        if (haploidize) {
                            sizes.add(popsampleposlist.get(i).size() - maxmissing[i] + 1);
                        } else {
                            sizes.add(2 * (popsampleposlist.get(i).size() - maxmissing[i]) + 1);
                        }
//                        System.out.println(sizes.get(i));
                        arraysize = arraysize.multiply(BigInteger.valueOf(sizes.get(i)));
                    }
                    dim = sizes.toArray(new Integer[0]);
                    SFS = new Integer[arraysize.intValue()];
                    Arrays.fill(SFS, 0);
//                    System.out.println(SFS.length);
//                    SFS = new int[2 * (pop1list.size() - maxmissing1) + 1][2 * (pop2list.size() - maxmissing2) + 1];
                    String headerout = null;
                    headerout = "m" + arrayofpop[0] + "\t" + "n" + arrayofpop[0];
                    for (int i = 1; i < arrayofpop.length; i++) {
                        headerout = headerout + "\t" + "m" + arrayofpop[i] + "\t" + "n" + arrayofpop[i];
                    }
                    writer.write("chromosome" + "\t" + "physpos" + "\t" + "genpos" + "\t" + headerout);
                    writer.newLine();
                }
            }
            writer.flush();
            writer.close();
        } catch (IOException ex) {
            Logger.getLogger(Main.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static final Pattern P_TAB = Pattern.compile("\t");
    private static final Pattern P_SPACE = Pattern.compile(" ");
    private static final Pattern P_COMMA = Pattern.compile(",");
    private static final boolean DEBUG = false;

    private static void reduceArray(List<Integer> list, Random r, int cap) {
        while (list.size() > cap) {
            int index = r.nextInt(list.size());
            list.remove(index);
        }
    }

    //extract 0 based indeces from a position in the one dimensional array
}
