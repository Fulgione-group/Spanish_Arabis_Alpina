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
//option handling needs to be modified!
public class Downsample_SFS {

    public static void main(String[] args) {
        System.out.println("Args: " + Arrays.toString(args));
//        downsampledSFS(new String[]{"--haploidize", "--vcf", "/home/btjeng/Data/testing/vcfShore_alpina810_2022-01-26.vcf.filteredQ30LD5UD100K.NoCaucasica2.SNPable.HW.mont.test.vcf", "--popfile", "/home/btjeng/Data/testing/alpina_all_samples_MOSTRECENT.poplabels", "--listofpop", "E6", "--ref", "mont", "--cutoff", "0.2"});
        downsampledSFS(new String[]{"--haploidize", "--vcf", "/home/btjeng/Data/testing/test_vcf_SFS.vcf", "--popfile", "/home/btjeng/Data/testing/1000Genomes_final_Canabria_3pxpclr.popfile", "--listofpop", "Group2,Group1,FRANCE", "--cutoff", "0.1"});
//        downsampledSFS(new String[]{"--vcf", "/home/btjeng/Data/testing/test_vcf_SFS.vcf", "--popfile", "/home/btjeng/Data/testing/1000Genomes_final_Canabria_3pxpclr.popfile", "--listofpop", "Group1,Group2,FRANCE", "--cutoff", "0.1"});
//        downsampledSFS(new String[]{"--haploidize","--vcf", "/home/btjeng/Data/testing/vcfShore_alpina810_2022-01-26.vcf.filteredQ30LD5UD100K.NoCaucasica2.SNPable.HW.mont.test.vcf", "--popfile", "/home/btjeng/Data/testing/alpina_all_samples_MOSTRECENT.poplabels", "--listofpop", "E3,E6", "--cutoff", "0.2"});
//                downsampledSFS("/home/btjeng/Data/testing/vcfShore_alpina810_2022-01-26.vcf.filteredQ30LD5UD100K.NoCaucasica2.SNPable.HW.mont.test.vcf", "/home/btjeng/Data/testing/alpina_all_samples_MOSTRECENT.poplabels", "E3", "E5", "mont", "/home/btjeng/Data/testing/test_out_SFS", Float.parseFloat("0.2"));
        // to build java byte code, all options are mandetory at the moment
//        downsampledSFS(args[0], args[1], args[2], args[3], args[4], args[5], Float.parseFloat(args[6]));
//        downsampledSFS(args);
    }

    // vcfpath1 = full path to the vcf
    // popfile = full path to popfile specifying a population ID for every sample ID separated by spaces. Each sample is writen on one line.
    // pop1 = population ID of population1
    // pop2 = population ID of population2
    // ref = sample ID of the outgroup
    // output = full path of the outputfile including its name
    private static void downsampledSFS(String[] inputList) {
        String inpath = null;
        String popfile = null;
        String listofpop = null;
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
            if (inputList[i].equals("--haploidize") | inputList[i].equals("-h")) {
                haploidize = true;
            }

        }
        if (outpath == null) {
            if (inpath.endsWith(".gz")) {
                outpath = inpath.substring(0, inpath.lastIndexOf(".")) + "." + listofpopout + ".SFS";
            } else {
                outpath = inpath + "." + listofpopout + ".SFS";
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
        String refsplit = null;
        String line;
        for (int i = 0; i < arrayofpop.length; i++) {
            popsamplelist.add(new ArrayList<String>());
            popsampleposlist.add(new ArrayList<Integer>());
            popsamplepostmp.add(new ArrayList<Integer>());
        }
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
            int helparray[] = {0, 2};
            while ((line = reader.readLine()) != null) {
                List<String> pop1tmp = new ArrayList<>();
                List<String> pop2tmp = new ArrayList<>();
                //Disregard header lines
                if (line.charAt(0) != '#') {
                    linearray = P_TAB.split(line);
                    if (linearray[4].contains(",")) {
//                    System.out.println(linearray[4]);
                        continue;
                    }
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
//                                System.out.println(linearray[popsampleposlist.get(i).get(j)].charAt(0));
                            } else {
                                popsamplepostmp.get(i).add(popsampleposlist.get(i).get(j));
//                                System.out.println(linearray[popsampleposlist.get(i).get(j)]);
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
//                            System.out.println(popsampleposlist.get(i).size() + "\t" + popsamplepostmp.get(i).size());
                            for (int j = 0; j < popsamplepostmp.get(i).size(); j++) {
                                genoarray = linearray[popsamplepostmp.get(i).get(j)];
//                                System.out.print(genoarray + "\t");
                                if (haploidize) {
                                    int randomOfTwoInts = r.nextBoolean() ? 0 : 2;
//                                    System.out.println(randomOfTwoInts);
                                    if (genoarray.charAt(randomOfTwoInts) != refsplit.charAt(0)) {
                                        acount[i]++;
//                                        System.out.println(genoarray.toString());
//                                        System.out.println(refsplit);
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
//                            System.out.println();
                        }
                        //add to two dimensional SFS
                        Integer index = mdarrayVal(dim, acount);
                        SFS[index]++;
//                        if (acount[0] == 0 & acount[1] == 0 & acount[2]==1) {
////                            System.out.println(Arrays.toString(dim));
//                            System.out.println(Arrays.toString(acount));
////                            System.out.println(popsamplepostmp.get(0).toString());
//                            System.out.println(index);
//                            System.out.println(Arrays.toString(mdarrayIndex(dim, index)));
//                        }

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
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(Main.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        Integer[] outarray = null;
        try ( BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath))))) {
            //output the SFS
            for (int i = 0; i < arrayofpop.length; i++) {
                writer.write(arrayofpop[i] + "\t");
            }
            writer.write("count");
            writer.newLine();
            for (int i = 0; i < SFS.length; i++) {
                outarray = mdarrayIndex(dim, i);
                Integer index_check = mdarrayVal(dim, outarray);
//                if(index_check == i){
//                    System.out.println(i + "\t" + index_check);
//                }
                for (int j = 0; j < outarray.length; j++) {
                    writer.write(outarray[j] + "\t");
                }
                writer.write("" + SFS[i]);
//                System.out.println(SFS[i]);
                writer.newLine();
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

    //get position from zero based index
    private static Integer mdarrayVal(Integer[] dim, Integer[] index) {
        Integer pos = -1;
        if (dim.length != index.length) {
            throw new IllegalArgumentException("dim length does not match indexing");
        } else {
            pos = index[0];
            for (int i = 1; i < dim.length; i++) {
                int dimfac = 1;
                for (int j = 0; j < i; j++) {
                    dimfac = dimfac * dim[j];
//                    System.out.println(dimfac);
                }
                pos = pos + index[i] * dimfac;
            }
        }
        return pos;
    }

    //extract 0 based indeces from a position in the one dimensional array
    private static Integer[] mdarrayIndex(Integer[] dim, Integer pos) {
        Integer[] index = new Integer[dim.length];
        Integer[] boundaries = new Integer[dim.length];
        Integer bound = 1;
        for (int i = 0; i < dim.length; i++) {
            bound = bound * dim[i];
            boundaries[i] = bound;
        }
//        System.out.println(Arrays.toString(boundaries));
        if (pos > bound) {
            throw new IllegalArgumentException("pos to big for set of dimension");
        } else {
//            System.out.println(boundaries[0]);
            index[0] = pos % boundaries[0];
            for (int i = 1; i < boundaries.length; i++) {
//                    index[i] = (pos / boundaries[(i-1)]) + ((pos % boundaries[(i-1)])/(pos % boundaries[(i-1)]))  - ((pos / boundaries[i]) - 1 + ((pos % boundaries[(i)])/(pos % boundaries[(i)]))) * boundaries[i];
                index[i] = (pos - (pos / boundaries[i]) * boundaries[i]) / boundaries[(i - 1)];
//                index[i] = (pos / boundaries[(i - 1)]) ;
//                System.out.println(((pos / boundaries[i]) - 1) * boundaries[i]);
            }
        }
//        System.out.println("Pos\t" + pos + "\tIndex\t" + Arrays.toString(index) + "\tboundaries\t" + Arrays.toString(boundaries));
        return index;
    }
}
