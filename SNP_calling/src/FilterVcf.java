/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Evolutionary_genetics;

import static Evolutionary_genetics.Main.flipZeroOne;
import static Evolutionary_genetics.Main.replaceCharAt;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author btjeng
 */
public class FilterVcf {

    public static void main(String[] args) {
        System.out.println("Args: " + Arrays.toString(args));
//        Filter_vcf(new String[]{"--vcf", "/home/btjeng/Data/MPIshore5/vcfShore_alpina810_2022-01-26.vcf.filteredQ30LD5UD100K_hcorrected_test.vcf", "--out", "/home/btjeng/Data/MPIshore5/vcfShore_alpina810_2022-01-26.vcf.filteredQ30LD5UD100K_hcorrected_test.nocaucasica.vcf", "--exclude", "/home/btjeng/Data/MPIshore5/exclude_MPIshore5.txt"});
//        Filter_vcf(new String[]{"--vcf", "/home/btjeng/Data/testing/GATK4.2_all_chr7.wheader.failed.vcf.gz", "--out", "/home/btjeng/Data/testing/failed_chr7.out", "--coverage", "/home/btjeng/Data/testing/GATK4.2_all_chrall2.GW_avcov.txt,-1", "--depth", "5,100000", "--min-sample-cov", "10", "--gzip", "--indels-exclude", "--quality", "30", "--biallelic"});
// use this to build byte code
        Filter_vcf(args);
    }

    private static void Filter_vcf(String[] inputList) {
        String inpath = null;
        String outpath = null;
        String outgroup = null;
        String excludeSampleFile = null;
        BufferedReader readerPop = null;
        BufferedWriter writer = null;
        boolean compress = false;
        boolean indels = false;
        boolean biallelic = false;
        boolean quality = false;
        boolean exclude_low_cov = false;
        boolean exclude_samples = false;
        boolean polarize = false;
        boolean checkavcov = false;
        int qthresh = 0;
        boolean depth = false;
        int dlthresh = 0;
        int duthresh = 0;
        float covfactor = 0;
        int minsamplecov = 0;
        int posCounter = 0;
        int retainPosCounter = 0;
        String meancovfile = null;
        Map<Integer, Float> mean_coverage = new HashMap<>();
        Map<String, Integer> genotype_pos = new HashMap<>();
        ArrayList<Integer> excludeList = new ArrayList<>();
        ArrayList<Integer> keepgenotypes = new ArrayList<Integer>();
        ArrayList<Integer> keeppos = new ArrayList<Integer>();
        for (int i = 0; i < inputList.length; i++) {
            // full input path of the vcf file
            System.out.println(i);
            if (inputList[i].equals("--vcf") | inputList[i].equals("-i")) {
                inpath = inputList[(i + 1)];
                i++;
            }
            //full path of the output vcf 
            if (inputList[i].equals("--out") | inputList[i].equals("-o")) {
                outpath = inputList[(i + 1)];
                i++;
            }
            // if specified indels will be excluded
            if (inputList[i].equals("--indels-exclude")) {
                indels= true;
            }
            // if specified only biallelic SNPs will be outputed along with reference calls
            if (inputList[i].equals("--biallelic")) {
                biallelic = true;
            }
            // if specified all genotypes with less quality then specified will be set to NA
            if (inputList[i].equals("--quality")) {
                quality = true;
                qthresh = Integer.parseInt(inputList[(i + 1)]);
                i++;
            }
            //everything below given minimum and given maximum coverage will be excluded
            //takes two arguments seperated by comma 
            if (inputList[i].equals("--depth")) {
                depth = true;
                dlthresh = Integer.parseInt(inputList[(i + 1)].split(",")[0]);
                duthresh = Integer.parseInt(inputList[(i + 1)].split(",")[1]);
                i++;
            }
            //if specified it will read average coverage per sample and will set every genotype with coverage greater than covfactor * average to NA
            //takes two arguments seperated by comma
            if (inputList[i].equals("--coverage")) {
                depth = true;
                meancovfile = inputList[(i + 1)].split(",")[0];
                covfactor = Float.parseFloat(inputList[(i + 1)].split(",")[1]);
                if (covfactor > 0) {
                    checkavcov = true;
                }
                i++;
            }
            if (inputList[i].equals("--exclude") | inputList[i].equals("-e")) {
                exclude_samples = true;
                excludeSampleFile = inputList[(i + 1)];
                i++;
            }
            //if specified it will filter genotypes based on average coverage from the coverage file
            //only works when --coverage is active
            if (inputList[i].equals("--min-sample-cov") & meancovfile != null) {
                exclude_low_cov = true;
                minsamplecov = Integer.parseInt(inputList[(i + 1)]);
                i++;
            }
            // if specified positions with heterozygout outgroup calls will be omitted and not outputted.
            // if the outgroup is derived in respect to the reference the SNPs will be polarized
            if (inputList[i].equals("--outgroup") | inputList[i].equals("-og")) {
                polarize = true;
                //needs to be the exact same as the sample name in the vcf
                outgroup = inputList[(i + 1)];
                i++;
            }
            // if specified output vcf will be gzipped
            if (inputList[i].equals("--gzip") | inputList[i].equals("-c")) {
                compress = true;
            }
        }
        if (outpath == null) {
            if (inpath.endsWith(".gz")) {
                outpath = inpath.substring(0, inpath.lastIndexOf(".")) + ".filtered.vcf";
            } else {
                outpath = inpath + ".filtered.vcf";
            }
        }
        if (inpath == null) {
            throw new IllegalArgumentException("No input specified");
        }
        try {
            // create input and output streams with regards to compression
            if (compress) {
                outpath = outpath + ".gz";
                writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(outpath)))));
            } else {
                writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outpath))));
            }
            InputStream inputStream = new BufferedInputStream(new FileInputStream(inpath));
            if (inpath.endsWith(".gz")) {
                inputStream = new GZIPInputStream(inputStream);
            }
            BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
            String[] linearray;
            String line;
            String[] genoline = null;
            //read mean cov for each sample
            //checks if meancov file was set
            if (meancovfile != null | excludeSampleFile != null) {
                InputStream inputStream2 = new BufferedInputStream(new FileInputStream(inpath));
                if (inpath.endsWith(".gz")) {
                    inputStream2 = new GZIPInputStream(inputStream2);
                }
                BufferedReader reader2 = new BufferedReader(new InputStreamReader(inputStream2));
                //go through the file until the header with sample IDs is found
                while ((line = reader2.readLine()) != null) {
                    if (line.startsWith("#CHROM")) {
                        linearray = P_TAB.split(line);
                        for (int i = 9; i < linearray.length; i++) {
                            genotype_pos.put(linearray[i], i);
                        }
                        genoline = linearray;
                        break;
                    }
                }
                if (meancovfile != null) {
                    InputStream inputStreamcovfile = new BufferedInputStream(new FileInputStream(meancovfile));
                    BufferedReader readercovfile = new BufferedReader(new InputStreamReader(inputStreamcovfile));
                    //skip first line (header)
                    line = readercovfile.readLine();
                    // tab separeted file with header where first column is the sample name as in the vcfheader and second column the mean coverage of the respective sample
                    while ((line = readercovfile.readLine()) != null) {
                        linearray = P_TAB.split(line);
                        mean_coverage.put(genotype_pos.get(linearray[0]), Float.parseFloat(linearray[1]));
                    }
                }
                if (excludeSampleFile != null) {
                    InputStream inputStreamExcludeFile = new BufferedInputStream(new FileInputStream(excludeSampleFile));
                    BufferedReader readerExcludeFile = new BufferedReader(new InputStreamReader(inputStreamExcludeFile));
                    //each genotype that needs to be excluded needs to be 
                    while ((line = readerExcludeFile.readLine()) != null) {
                        excludeList.add(genotype_pos.get(line));
                    }
                }
                //indentify low coverage samples and print if samples are included or excluded
                //if out group is excluded  set polarize to false
                if (exclude_low_cov || exclude_samples) {
                    for (int i = 0; i < 9; i++) {
                        keeppos.add(i);
                    }
                    for (int i = 9; i < genoline.length; i++) {
//                        System.out.println(mean_coverage.get(i));
                        if (exclude_low_cov && exclude_samples) {
                            if (mean_coverage.get(i) > minsamplecov && !excludeList.contains(i)) {
                                keepgenotypes.add(i);
                                keeppos.add(i);
                                System.out.println("Include genotype:" + "\t" + genoline[i] + "\tCoverage:\t" + mean_coverage.get(i) + "\tfor threshold\t" + minsamplecov);
                            } else {
                                if (genoline[i] == outgroup) {
                                    System.out.println("Warning outgroup has too low mean coverage. VCF will not be polarized!. Coverage:\t");
                                    polarize = false;
                                }
                                System.out.println("Exclude genotype, coverage or --exclude " + excludeSampleFile + ":\t" + genoline[i] + "\tCoverage:\t" + mean_coverage.get(i) + "\tfor threshold\t" + minsamplecov);
                            }
                        } else if (exclude_low_cov) {
                            if (mean_coverage.get(i) > minsamplecov) {
                                keepgenotypes.add(i);
                                keeppos.add(i);
                                System.out.println("Include genotype:" + "\t" + genoline[i] + "\tCoverage:\t" + mean_coverage.get(i) + "\tfor threshold\t" + minsamplecov);
                            } else {
                                if (genoline[i] == outgroup) {
                                    System.out.println("Warning outgroup has too low mean coverage. VCF will not be polarized!. Coverage:\t");
                                    polarize = false;
                                }
                                System.out.println("Exclude genotype:" + "\t" + genoline[i] + "\tCoverage:\t" + mean_coverage.get(i) + "\tfor threshold\t" + minsamplecov);
                            }
                        } else {
                            if (!excludeList.contains(i)) {
                                keepgenotypes.add(i);
                                keeppos.add(i);
                                System.out.println("Include genotype:" + "\t" + genoline[i]);
                            } else {
                                if (genoline[i] == outgroup) {
                                    System.out.println("Warning outgroup has too low mean coverage. VCF will not be polarized!. Coverage:\t");
                                    polarize = false;
                                }
                                System.out.println("Exclude genotype --exclude " + excludeSampleFile + ":\t" + genoline[i]);
                            }
                        }
                    }
                }
            }
            String[] format;
            int nonmissing;
            int linecounter = 0;
            String[] genotype;
            int outgrouppos = 0;
            String[] outgroup_geno = null;
            boolean calledsite = false;
            while ((line = reader.readLine()) != null) {
                linecounter++;
                //omit header lines
                if (!line.startsWith("#")) {
                    System.out.println(line);
                    posCounter++;
                    linearray = P_TAB.split(line);
                    format = linearray[8].split(":");
                    int GQpos = -1;
                    int DPpos = -1;
//                    System.out.println(linearray[8]);
//                    System.out.println(Arrays.asList(format));
                    calledsite = false;
                    for (int i = 0; i < format.length; i++) {
                        if (format[i].equals("GQ") | format[i].equals("RGQ")) {
                            GQpos = i;
                        }
                        if (format[i].equals("DP")) {
                            DPpos = i;
                        }
                    }
//                    System.out.println(GQpos +"\t" + DPpos);
                    if (indels) {
                        if (linearray[3].length() > 1 || linearray[4].length() > 1 || linearray[4].equals("-") || linearray[4].equals("*")) {
                            continue;
                        }
                    }
                    //if in alt is a comma to sep the to alleles the position will be omitted
                    if (biallelic) {
                        if (linearray[4].contains(",")) {
                            continue;
                        }
                    }
                    //if outgroup is specified
                    if (polarize) {
                        boolean flipref = false;
                        outgroup_geno = linearray[outgrouppos].split(":");
                        boolean genEqualsPipe = outgroup_geno[0].equals(".|.");
                        boolean genEqualsSlash = outgroup_geno[0].equals("./.");
                        boolean genEqualsHet = (outgroup_geno[0].charAt(0) != outgroup_geno[0].charAt(2));
                        boolean genLength = (outgroup_geno.length == 1);
                        //check if genotype is not NA in the first place
                        if (genEqualsPipe || genEqualsSlash || genEqualsHet || genLength) {
                            continue;
                        } else {
                            int parsedGQpos = Integer.parseInt(outgroup_geno[GQpos]);
                            int parsedDPpos = Integer.parseInt(outgroup_geno[DPpos]);
                            if (quality && depth) {
                                if (checkavcov) {
                                    //Set genotypes to NA if they do not meet every of these metrics. 
                                    //if coverage and depth are specified it will take the most conservative threshold
                                    if (parsedGQpos < qthresh || parsedDPpos < dlthresh || parsedDPpos > duthresh || parsedDPpos > covfactor * mean_coverage.get(outgrouppos)) {
                                        continue;
                                    } else {
                                        // if the outgroup is derived in regard to the reference sequence the reference and the alternative allele need to be flipped 
                                        if (outgroup_geno[0].charAt(0) == '1') {
                                            flipref = true;
                                        }
                                    }
                                } else {
                                    if (parsedGQpos < qthresh || parsedDPpos < dlthresh || parsedDPpos > duthresh) {
                                        continue;
                                    } else {
                                        if (outgroup_geno[0].charAt(0) == '1') {
                                            flipref = true;
                                        }
                                    }
                                }
                            } else if (depth) {
                                if (checkavcov) {
                                    if (parsedDPpos < dlthresh || parsedDPpos > duthresh || parsedDPpos > covfactor * mean_coverage.get(outgrouppos)) {
                                        continue;
                                    } else {
                                        if (outgroup_geno[0].charAt(0) == '1') {
                                            flipref = true;
                                        }
                                    }
                                } else {
                                    if (parsedDPpos < dlthresh || parsedDPpos > duthresh) {
                                        continue;
                                    } else {
                                        if (outgroup_geno[0].charAt(0) == '1') {
                                            flipref = true;
                                        }
                                    }
                                }
                            } else if (quality) {
                                if (parsedGQpos < qthresh) {
                                    continue;
                                } else {
                                    if (outgroup_geno[0].charAt(0) == '1') {
                                        flipref = true;
                                    }
                                }
                            } else {
                                if (outgroup_geno[0].charAt(0) == '1') {
                                    flipref = true;
                                }
                            }
                        }
                        // flipp the annotation of the reference and alt allele in that position
                        if (flipref) {
                            String ref = linearray[3];
                            String alt = linearray[4];
                            linearray[3] = alt;
                            linearray[4] = ref;
                        }
                        if (quality & depth) {
                            for (int i = 9; i < linearray.length; i++) {
                                genotype = linearray[i].split(":");
                                if (!checkavcov) {
                                    genEqualsSlash = genotype[0].equals("./.");
                                    genEqualsPipe = genotype[0].equals(".|.");
                                    genLength = (genotype.length == 1);
                                    if (!genEqualsSlash && !genEqualsPipe) {
                                        if (genLength) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            int parsedGQpos = Integer.parseInt(genotype[GQpos]);
                                            int parsedDPpos = Integer.parseInt(genotype[DPpos]);
                                            if (parsedGQpos < qthresh || parsedDPpos < dlthresh || parsedDPpos > duthresh) {
                                                linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                                linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                            } else {
                                                calledsite = true;
                                                if (flipref) {
                                                    linearray[i] = flipZeroOne(linearray[i], 0);
                                                    linearray[i] = flipZeroOne(linearray[i], 2);
                                                }
                                            }
                                        }
                                    }
                                } else {
//                                System.out.println("i: " + i + ", line: " + linecounter);
                                    genEqualsSlash = genotype[0].equals("./.");
//                                System.out.println("genotype: " + Arrays.asList(genotype) + ", GQpos: " + GQpos + "\t" + linearray[1]);
                                    genEqualsPipe = genotype[0].equals(".|.");
                                    genLength = (genotype.length == 1);
                                    if (!genEqualsSlash && !genEqualsPipe) {
                                        if (genLength) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            int parsedGQpos = Integer.parseInt(genotype[GQpos]);
                                            int parsedDPpos = Integer.parseInt(genotype[DPpos]);
//                                    System.out.println("genEquals: " + genEqualsSlash + ", genQualsPipe; " + genEqualsPipe + ", parsedGQpos: " + parsedGQpos + ", parsedDPpos: " + parsedDPpos);
                                            //System.out.println(Arrays.asList(genotype) + "\t" +genotype[GQpos] +"\t" +genotype[DPpos]+ "\t" +mean_coverage.get(i));
                                            if (parsedGQpos < qthresh || parsedDPpos < dlthresh || parsedDPpos > duthresh || parsedDPpos > covfactor * mean_coverage.get(i)) {
                                                linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                                linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                            } else {
                                                calledsite = true;
                                                if (flipref) {
                                                    linearray[i] = flipZeroOne(linearray[i], 0);
                                                    linearray[i] = flipZeroOne(linearray[i], 2);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        } else if (quality) {
                            for (int i = 9; i < linearray.length; i++) {
                                genotype = linearray[i].split(":");
//                                System.out.println("i: " + i + ", line: " + linecounter);
                                genEqualsSlash = genotype[0].equals("./.");
//                                System.out.println("genotype: " + Arrays.asList(genotype) + ", GQpos: " + GQpos + "\t" + linearray[1]);
                                genEqualsPipe = genotype[0].equals(".|.");
                                genLength = (genotype.length == 1);
                                if (!genEqualsSlash && !genEqualsPipe) {
                                    if (genLength) {
                                        linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                        linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                    } else {
                                        int parsedGQpos = Integer.parseInt(genotype[GQpos]);
//                                    System.out.println("genEquals: " + genEqualsSlash + ", genQualsPipe; " + genEqualsPipe + ", parsedGQpos: " + parsedGQpos + ", parsedDPpos: " + parsedDPpos);
                                        //System.out.println(Arrays.asList(genotype) + "\t" +genotype[GQpos] +"\t" +genotype[DPpos]+ "\t" +mean_coverage.get(i));
                                        if (parsedGQpos < qthresh) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            // if all other genotypes passed QC but outgroup is derived flip the calls
                                            calledsite = true;
                                            if (flipref) {
                                                linearray[i] = flipZeroOne(linearray[i], 0);
                                                linearray[i] = flipZeroOne(linearray[i], 2);
                                            }
                                        }
                                    }
                                }
                            }
                        } else if (depth) {
                            for (int i = 9; i < linearray.length; i++) {
                                genotype = linearray[i].split(":");
                                if (!checkavcov) {
//                                System.out.println("i: " + i + ", line: " + linecounter);
                                    genEqualsSlash = genotype[0].equals("./.");
//                                System.out.println("genotype: " + Arrays.asList(genotype) + ", GQpos: " + GQpos + "\t" + linearray[1]);
                                    genEqualsPipe = genotype[0].equals(".|.");
                                    genLength = (genotype.length == 1);
                                    if (!genEqualsSlash && !genEqualsPipe) {
                                        if (genLength) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            int parsedDPpos = Integer.parseInt(genotype[DPpos]);
//                                    System.out.println("genEquals: " + genEqualsSlash + ", genQualsPipe; " + genEqualsPipe + ", parsedGQpos: " + parsedGQpos + ", parsedDPpos: " + parsedDPpos);
                                            //System.out.println(Arrays.asList(genotype) + "\t" +genotype[GQpos] +"\t" +genotype[DPpos]+ "\t" +mean_coverage.get(i));
                                            if (parsedDPpos < dlthresh || parsedDPpos > duthresh) {
                                                linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                                linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                            } else {
                                                calledsite = true;
                                                if (flipref) {
                                                    linearray[i] = flipZeroOne(linearray[i], 0);
                                                    linearray[i] = flipZeroOne(linearray[i], 2);
                                                }
                                            }
                                        }
                                    }
                                } else {
//                                System.out.println("i: " + i + ", line: " + linecounter);
                                    genEqualsSlash = genotype[0].equals("./.");
//                                System.out.println("genotype: " + Arrays.asList(genotype) + ", GQpos: " + GQpos + "\t" + linearray[1]);
                                    genEqualsPipe = genotype[0].equals(".|.");
                                    genLength = (genotype.length == 1);
                                    if (!genEqualsSlash && !genEqualsPipe) {
                                        if (genLength) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            int parsedDPpos = Integer.parseInt(genotype[DPpos]);
//                                    System.out.println("genEquals: " + genEqualsSlash + ", genQualsPipe; " + genEqualsPipe + ", parsedGQpos: " + parsedGQpos + ", parsedDPpos: " + parsedDPpos);
                                            //System.out.println(Arrays.asList(genotype) + "\t" +genotype[GQpos] +"\t" +genotype[DPpos]+ "\t" +mean_coverage.get(i));
                                            if (parsedDPpos < dlthresh || parsedDPpos > duthresh || parsedDPpos > covfactor * mean_coverage.get(i)) {
                                                linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                                linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                            } else {
                                                calledsite = true;
                                                if (flipref) {
                                                    linearray[i] = flipZeroOne(linearray[i], 0);
                                                    linearray[i] = flipZeroOne(linearray[i], 2);
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                        }
                        // no polarization                            
                    } else {
                        if (quality & depth) {
                            for (int i = 9; i < linearray.length; i++) {
                                genotype = linearray[i].split(":");
//                            System.out.println(Arrays.asList(genotype));
//                            System.out.println(GQpos);
//                            System.out.println(genotype[GQpos]);
//                            System.out.println(i);
                                if (!checkavcov) {
//                                System.out.println("i: " + i + ", line: " + linecounter);
                                    boolean genEqualsSlash = genotype[0].equals("./.");
//                                System.out.println("genotype: " + Arrays.asList(genotype) + ", GQpos: " + GQpos + "\t" + linearray[1]);
                                    boolean genEqualsPipe = genotype[0].equals(".|.");
                                    boolean genLength = (genotype.length == 1);
                                    if (!genEqualsSlash && !genEqualsPipe) {
                                        if (genLength) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            int parsedGQpos = Integer.parseInt(genotype[GQpos]);
                                            int parsedDPpos = Integer.parseInt(genotype[DPpos]);
//                                    System.out.println("genEquals: " + genEqualsSlash + ", genQualsPipe; " + genEqualsPipe + ", parsedGQpos: " + parsedGQpos + ", parsedDPpos: " + parsedDPpos);
                                            //System.out.println(Arrays.asList(genotype) + "\t" +genotype[GQpos] +"\t" +genotype[DPpos]+ "\t" +mean_coverage.get(i));
                                            if (parsedGQpos < qthresh || parsedDPpos < dlthresh || parsedDPpos > duthresh) {
                                                linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                                linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                            } else {
                                                calledsite = true;
                                            }
                                        }
                                    }
                                } else {
//                                System.out.println("i: " + i + ", line: " + linecounter);
                                    boolean genEqualsSlash = genotype[0].equals("./.");
//                                System.out.println("genotype: " + Arrays.asList(genotype) + ", GQpos: " + GQpos + "\t" + linearray[1]);
                                    boolean genEqualsPipe = genotype[0].equals(".|.");
                                    boolean genLength = (genotype.length == 1);
                                    if (!genEqualsSlash && !genEqualsPipe) {
                                        if (genLength) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            int parsedGQpos = Integer.parseInt(genotype[GQpos]);
                                            int parsedDPpos = Integer.parseInt(genotype[DPpos]);
//                                    System.out.println("genEquals: " + genEqualsSlash + ", genQualsPipe; " + genEqualsPipe + ", parsedGQpos: " + parsedGQpos + ", parsedDPpos: " + parsedDPpos);
                                            //System.out.println(Arrays.asList(genotype) + "\t" +genotype[GQpos] +"\t" +genotype[DPpos]+ "\t" +mean_coverage.get(i));
                                            if (parsedGQpos < qthresh || parsedDPpos < dlthresh || parsedDPpos > duthresh || parsedDPpos > covfactor * mean_coverage.get(i)) {
                                                linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                                linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                            } else {
                                                calledsite = true;
                                            }
                                        }
                                    }
                                }
                            }
                        } else if (quality) {
                            for (int i = 9; i < linearray.length; i++) {
                                genotype = linearray[i].split(":");
//                                System.out.println("i: " + i + ", line: " + linecounter);
                                boolean genEqualsSlash = genotype[0].equals("./.");
//                                System.out.println("genotype: " + Arrays.asList(genotype) + ", GQpos: " + GQpos + "\t" + linearray[1]);
                                boolean genEqualsPipe = genotype[0].equals(".|.");
                                boolean genLength = (genotype.length == 1);
                                if (!genEqualsSlash && !genEqualsPipe) {
                                    if (genLength) {
                                        linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                        linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                    } else {
                                        int parsedGQpos = Integer.parseInt(genotype[GQpos]);
//                                    System.out.println("genEquals: " + genEqualsSlash + ", genQualsPipe; " + genEqualsPipe + ", parsedGQpos: " + parsedGQpos + ", parsedDPpos: " + parsedDPpos);
                                        //System.out.println(Arrays.asList(genotype) + "\t" +genotype[GQpos] +"\t" +genotype[DPpos]+ "\t" +mean_coverage.get(i));
                                        if (parsedGQpos < qthresh) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            calledsite = true;
                                        }
                                    }
                                }
                            }
                        } else if (depth) {
                            for (int i = 9; i < linearray.length; i++) {
                                genotype = linearray[i].split(":");
                                if (!checkavcov) {
//                                System.out.println("i: " + i + ", line: " + linecounter);
                                    boolean genEqualsSlash = genotype[0].equals("./.");
//                                System.out.println("genotype: " + Arrays.asList(genotype) + ", GQpos: " + GQpos + "\t" + linearray[1]);
                                    boolean genEqualsPipe = genotype[0].equals(".|.");
                                    boolean genLength = (genotype.length == 1);
                                    if (!genEqualsSlash && !genEqualsPipe) {
                                        if (genLength) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            int parsedDPpos = Integer.parseInt(genotype[DPpos]);
//                                    System.out.println("genEquals: " + genEqualsSlash + ", genQualsPipe; " + genEqualsPipe + ", parsedGQpos: " + parsedGQpos + ", parsedDPpos: " + parsedDPpos);
                                            //System.out.println(Arrays.asList(genotype) + "\t" +genotype[GQpos] +"\t" +genotype[DPpos]+ "\t" +mean_coverage.get(i));
                                            if (parsedDPpos < dlthresh || parsedDPpos > duthresh) {
                                                linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                                linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                            } else {
                                                calledsite = true;
                                            }
                                        }
                                    }
                                } else {
//                                System.out.println("i: " + i + ", line: " + linecounter);
                                    boolean genEqualsSlash = genotype[0].equals("./.");
//                                System.out.println("genotype: " + Arrays.asList(genotype) + ", GQpos: " + GQpos + "\t" + linearray[1]);
                                    boolean genEqualsPipe = genotype[0].equals(".|.");
                                    boolean genLength = (genotype.length == 1);
                                    if (!genEqualsSlash && !genEqualsPipe) {
                                        if (genLength) {
                                            linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                            linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                        } else {
                                            int parsedDPpos = Integer.parseInt(genotype[DPpos]);
//                                    System.out.println("genEquals: " + genEqualsSlash + ", genQualsPipe; " + genEqualsPipe + ", parsedGQpos: " + parsedGQpos + ", parsedDPpos: " + parsedDPpos);
                                            //System.out.println(Arrays.asList(genotype) + "\t" +genotype[GQpos] +"\t" +genotype[DPpos]+ "\t" +mean_coverage.get(i));
                                            if (parsedDPpos < dlthresh || parsedDPpos > duthresh || parsedDPpos > covfactor * mean_coverage.get(i)) {
                                                linearray[i] = replaceCharAt(linearray[i], 0, '.');
                                                linearray[i] = replaceCharAt(linearray[i], 2, '.');
                                            } else {
                                                calledsite = true;
                                            }
                                        }
                                    }
                                }
                            }

                        } else {
                            calledsite = true;
                        }
                    }
                    // if non of the sites passed QC exlcude that position from the file
                    if (!calledsite) {
                        continue;
                    }
                    retainPosCounter++;
                    String tmpline = linearray[0];
                    if (exclude_low_cov || exclude_samples) {
                        for (int i = 1; i < keeppos.size(); i++) {
                            tmpline = tmpline + "\t" + linearray[keeppos.get(i)];
                        }
                    } else {
                        for (int i = 1; i < linearray.length; i++) {
                            tmpline = tmpline + "\t" + linearray[i];
                        }
                    }
                    writer.write(tmpline);
                    writer.newLine();
                } else if (line.startsWith("#CHROM")) {
                    linearray = P_TAB.split(line);
                    if (polarize) {
                        for (int i = 9; i < linearray.length; i++) {
                            if (outgroup.equals(linearray[i])) {
                                outgrouppos = i;
                            }
                        }
                    }
                    if (exclude_low_cov || exclude_samples) {
                        String tmpline = linearray[0];
                        for (int i = 1; i < keeppos.size(); i++) {
                            tmpline = tmpline + "\t" + linearray[keeppos.get(i)];
                        }
                        System.out.println(keeppos.size());
                        writer.write(tmpline);
                        writer.newLine();
                    } else {
                        writer.write(line);
                        writer.newLine();
                    }
                } else {
                    writer.write(line);
                    writer.newLine();
                }
            }
            //flush and close the output stream to avoid broken file end
            writer.flush();
            writer.close();
            System.out.println("total_positions_before:\t" + posCounter + "\ttotal_positions_after:\t" + retainPosCounter);
        } catch (IOException ex) {
            Logger.getLogger(Main.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
    }

    //replace character in a String used to set genotypes that fail QC to NA
    public static String replaceCharAt(String s, int pos, char c) {
        return s.substring(0, pos) + c + s.substring(pos + 1);
    }

    //flip zeros to ones and ones to zeros.
    public static String flipZeroOne(String s, int pos) {
        return s.substring(0, pos) + (1 - Integer.valueOf(String.valueOf(s.charAt(pos)))) + s.substring(pos + 1);
    }
    //Split a string into an Array fast
    private static final Pattern P_TAB = Pattern.compile("\t");

}
