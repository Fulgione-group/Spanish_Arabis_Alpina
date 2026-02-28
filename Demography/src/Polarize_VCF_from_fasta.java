package Evolutionary_genetics;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author btjeng
 */
public class Polarize_VCF_from_fasta {

    public static void main(String[] args) throws IOException {
//        System.out.println("Args: " + Arrays.toString(args));
//        Polarize(new String[]{"--vcf", "/home/btjeng/Data/testing/test.vcf", "--fasta", "/home/btjeng/Data/testing/all_withPatch_merged_noID_sorted_alpina_syntenic_alpinaRef_mont_zeroBasedStart_oneBasedEnd_chrall_sorted_mlines_noGaps_uppercases.fasta"});
        Polarize(args);
    }

    private static void Polarize(String[] inputList) throws IOException {
        String vcfpath = null;
        String fastapath = null;
        String outpath = null;
        String line = null;
        String fasta_contig = null;
        int fasta_contig_counter = -1;
        int vcf_contig_counter = -1;
        int contig_pos = -1;
        int linecounter = 0;
        boolean stdout = false;
        ArrayList<ArrayList<Character>> anc_states = new ArrayList<>();
        for (int i = 0; i < inputList.length; i++) {
            // full input path of the vcf file
            if (inputList[i].equals("--vcf") | inputList[i].equals("-i")) {
                vcfpath = inputList[(i + 1)];
                i++;
            }//contig ofs fasta need to have the same order as
            if (inputList[i].equals("--fasta") | inputList[i].equals("-i")) {
                fastapath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--out") | inputList[i].equals("-i")) {
                outpath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--stdout") | inputList[i].equals("-i")) {
                stdout = true;
            }
        }
        InputStream vcf_inputStream = new BufferedInputStream(new FileInputStream(vcfpath));
        if (vcfpath.endsWith(".gz")) {
            vcf_inputStream = new GZIPInputStream(vcf_inputStream);
            if (outpath == null) {
                outpath = vcfpath.substring(0, vcfpath.lastIndexOf(".")) + ".polarized.vcf";
            }
        } else {
            if (outpath == null) {
                outpath = vcfpath + ".polarized.vcf";
            }
        }
//        System.out.println(outpath);
        BufferedReader vcf_reader = new BufferedReader(new InputStreamReader(vcf_inputStream));
        InputStream fasta_inputStream = new BufferedInputStream(new FileInputStream(fastapath));
        if (fastapath.endsWith(".gz")) {
            fasta_inputStream = new GZIPInputStream(fasta_inputStream);
        }
        BufferedReader fasta_reader = new BufferedReader(new InputStreamReader(fasta_inputStream));
        while ((line = fasta_reader.readLine()) != null) {
            if (line.startsWith(">")) {
                fasta_contig = line;
                fasta_contig_counter++;
                ArrayList<Character> contig_array = new ArrayList<>();
                anc_states.add(contig_array);
            } else {
                for (int i = 0; i < line.length(); i++) {
                    anc_states.get(fasta_contig_counter).add(line.charAt(i));
                }
            }
        }
        String contig = "start";
        String[] linearray = null;
        String tmpline = null;
        BufferedWriter writerVCF = null;
        if (!stdout) {
            writerVCF = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath))));
        }
        int pos = 0;
        while ((line = vcf_reader.readLine()) != null) {
            linecounter++;
            if (!line.startsWith("#")) {
                linearray = P_TAB.split(line);
                if (!contig.equals(linearray[0])) {
                    contig = linearray[0];
                    vcf_contig_counter++;
//                    System.out.println(contig);
                }
                pos = (Integer.valueOf(linearray[1]) - 1);
//                boolean isN = (anc_states.get(vcf_contig_counter).get((pos - 1)) == 'N');
//                System.out.println(linearray[3].charAt(0));
//                System.out.println(anc_states.get(vcf_contig_counter).get((pos - 1)));
                boolean isMultiAllelic = linearray[4].contains(",");
                if(isMultiAllelic){
                    continue;
                }
                boolean isInvariantRef = (linearray[4].charAt(0)=='.');
                boolean isN = (anc_states.get(vcf_contig_counter).get(pos)=='N');
                boolean isAnc = (anc_states.get(vcf_contig_counter).get(pos) == linearray[3].charAt(0));
                boolean isDir = (anc_states.get(vcf_contig_counter).get(pos) == linearray[4].charAt(0));
                if (isDir) {
                    String Anc = linearray[4];
                    String Dir = linearray[3];
                    linearray[3] = Anc;
                    linearray[4] = Dir;
                    linearray[6] = "anc_flip";
                    for (int i = 9; i < linearray.length; i++) {
                        boolean isRef1 = (linearray[i].charAt(0) == '0');
                        boolean isAlt1 = (linearray[i].charAt(0) == '1');
                        boolean isRef2 = (linearray[i].charAt(2) == '0');
                        boolean isAlt2 = (linearray[i].charAt(2) == '1');
                        if (isRef1) {
                            linearray[i] = replaceCharAt(linearray[i], 0, '1');
                        }
                        if (isAlt1) {
                            linearray[i] = replaceCharAt(linearray[i], 0, '0');
                        }
                        if (isRef2) {
                            linearray[i] = replaceCharAt(linearray[i], 2, '1');
                        }
                        if (isAlt2) {
                            linearray[i] = replaceCharAt(linearray[i], 2, '0');
                        }
                    }
                    tmpline = linearray[0];
                    for (int i = 1; i < linearray.length; i++) {
                        tmpline = tmpline + "\t" + linearray[i];
                    }
                } else if(isInvariantRef && !isAnc && !isN){
                    String Anc = anc_states.get(vcf_contig_counter).get(pos).toString();
                    String Dir = linearray[3];
                    linearray[3] = Anc;
                    linearray[4] = Dir;
                    linearray[6] = "anc_flip";
                    for (int i = 9; i < linearray.length; i++) {
                        boolean isRef1 = (linearray[i].charAt(0) == '0');
                        boolean isAlt1 = (linearray[i].charAt(0) == '1');
                        boolean isRef2 = (linearray[i].charAt(2) == '0');
                        boolean isAlt2 = (linearray[i].charAt(2) == '1');
                        if (isRef1) {
                            linearray[i] = replaceCharAt(linearray[i], 0, '1');
                        }
                        if (isAlt1) {
                            linearray[i] = replaceCharAt(linearray[i], 0, '0');
                        }
                        if (isRef2) {
                            linearray[i] = replaceCharAt(linearray[i], 2, '1');
                        }
                        if (isAlt2) {
                            linearray[i] = replaceCharAt(linearray[i], 2, '0');
                        }
                    }
                    tmpline = linearray[0];
                    for (int i = 1; i < linearray.length; i++) {
                        tmpline = tmpline + "\t" + linearray[i];
                    }
                }
                if (isDir | (isInvariantRef && !isAnc && !isN)) {
                    if (!stdout) {
                        writerVCF.write(tmpline);
                        writerVCF.newLine();
                    } else {
                        System.out.println(tmpline);
                    }
                } else if (isAnc) {
                    if (!stdout) {
                        writerVCF.write(line);
                        writerVCF.newLine();
                    } else {
                        System.out.println(line);
                    }
                }
            } else {
                if (!stdout) {
                    writerVCF.write(line);
                    writerVCF.newLine();
                } else {
                    System.out.println(line);
                }
                if (linecounter == 1) {
                    if (!stdout) {
                        writerVCF.write("##FILTER=<ID=anc_flip,Description=\"this vcf is polarized. Anc state is changed contains anc_flip\">");
                        writerVCF.newLine();
                    } else {
                        System.out.println("##FILTER=<ID=anc_flip,Description=\"this vcf is polarized. Anc state is changed contains anc_flip\">");
                    }
                }
            }
        }
        writerVCF.flush();
        writerVCF.close();
    }

    public static String replaceCharAt(String s, int pos, char c) {
        return s.substring(0, pos) + c + s.substring(pos + 1);
    }
    private static final Pattern P_TAB = Pattern.compile("\t");

}
