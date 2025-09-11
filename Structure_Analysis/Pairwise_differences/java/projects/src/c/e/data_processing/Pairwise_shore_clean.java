package c.e.data_processing;

import java.io.File;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.*;
import java.io.*;
import java.lang.reflect.Constructor;

public class Pairwise_shore_clean {
	
	public Pairwise_shore_clean() {}
	
	public void setMatrixFile(String matrixName, String accessionString, String results, String maskName){
		
		try {
                        int accessionIndex = Integer.parseInt(accessionString) + 3;
                        // int accession = Integer.parseInt(accessionString);

			System.out.println("src/c/e/data_processing/Pairwise_shore_clean.java");			
                        System.out.println("Running on sample: " + accessionIndex);





			////
			//	Load the specific mask
			////
			int numChr = 8;
			BitSet[] mask = new BitSet[numChr];
			for (int c=0; c<numChr; c++) {
				mask[c]= new BitSet();
			}
			File fileMask = new File(maskName);
			Scanner scannerMask = new Scanner(fileMask);
		        while ( scannerMask.hasNextLine() ) {
                        	String snp = scannerMask.nextLine();
                        	String[] splitSnp = snp.split("\t");
				// if (Integer.parseInt(splitSnp[2]) > Integer.parseInt(splitSnp[1]) && Integer.parseInt(splitSnp[3]) == 1 ) {
				 if (Integer.parseInt(splitSnp[2]) > Integer.parseInt(splitSnp[1])) {
					int chrMask = Character.getNumericValue(splitSnp[0].charAt(3)) - 1;
					int min = Integer.parseInt(splitSnp[1]);
					int max = Integer.parseInt(splitSnp[2]);
					for (int i=min; i<=max; i++) {
						mask[chrMask].set(i);
					}
				}
			}
			
			int cardin = 0;
			for (int c=0; c<numChr; c++) {
			//	// mask[c].and(maskRepeat[c]);
				cardin = cardin + mask[c].cardinality();
			}
			System.out.println("Cardinality: " + cardin);





	    		// Begin with the big matrix
	    		File matrix = new File(matrixName);
	    		Scanner scannerMatrix = new Scanner(matrix);
	    		
	    		// Just get ID order in the matrix from first line
	    		// And count iberians [Spanish + Portuguese]
			String idsUnsplit = scannerMatrix.nextLine();
	    		String[] splitIDs = idsUnsplit.split("\t");	
			


			// Build IndexesToGet, array of interesting indexes 
			
			int allPlants = splitIDs.length - 3;
			int howMany = allPlants;
			System.out.println("howMany: "+ howMany);
			
			int[] IndexesToGet = new int[howMany];

			System.out.println("Focus on: " + accessionIndex + " " + splitIDs[accessionIndex]);
			for (int i=3; i<splitIDs.length; i++) {	
				IndexesToGet[i-3] = i;
			}
		
	
			//
			// Now the real game
			// 

       		 	// Begin to calculate pairwise differences
        	
			int[] differences = new int[IndexesToGet.length];
			int[] length = new int[IndexesToGet.length];
			int chr =0;
			
		    	matrix = new File(matrixName);
		    	scannerMatrix = new Scanner(matrix);
		    	scannerMatrix.nextLine();
	    	
			while ( scannerMatrix.hasNextLine() ) {
				String snp = scannerMatrix.nextLine();
				String[] splitSnp = snp.split("\t");
				// System.out.println(splitSnp.length);
				
				// Just print where we are
        			if (chr != Integer.parseInt(splitSnp[0])) {
        				System.out.println("Chr: " + Integer.parseInt(splitSnp[0]));
        			}
        			chr = Integer.parseInt(splitSnp[0]);
        			int pos = Integer.parseInt(splitSnp[1]);
	       		 	//
	        		
       		 		if (!mask[chr-1].get(pos)) {
    		        		char base1 = splitSnp[accessionIndex].charAt(0);
	    	        		
    			       		for (int acc2=0; acc2<IndexesToGet.length; acc2++) {
    			       			char base2 = splitSnp[IndexesToGet[acc2]].charAt(0);
    			       			
    			       			if ( (base1 != 'N') && (base2 != 'N') ) {
    			       				length[acc2] = length[acc2] + 1;
    		        				if (base1 != base2) {
    			       					differences[acc2] = differences[acc2] + 1;
    			       				}
    			       			}
    		        		}
       	 			} 		
		        }
			// Calculate pi
			
			double[] pairwDiff = new double[IndexesToGet.length];
			for (int q=0; q<IndexesToGet.length; q++) {
				pairwDiff[q] = (double)(differences[q])/(double)(length[q]);
			}
			// System.out.println("Check: len: " + length[indexHLUI] + " diff: " + differences[indexHLUI] + " pi: " + pairwDiff[indexHLUI]);
			// Done
			
			
			
			// Write out
			
			Writer writer = new FileWriter(results + "pi" + accessionIndex + ".txt");
			PrintWriter out = new PrintWriter(writer);
			for (int q=0; q<IndexesToGet.length; q++) {
				out.print(pairwDiff[q]);
				out.print("\t");
			}
			out.print("\n");
			out.close();



		
                        Writer writer2 = new FileWriter(results + "len" + accessionIndex + ".txt");
                        PrintWriter out2 = new PrintWriter(writer2);
                        for (int q=0; q<IndexesToGet.length; q++) {
                                out2.print(length[q]);
                                out2.print("\t");
                        }
                        out2.print("\n");
                        out2.close();




                        Writer writer3 = new FileWriter(results + "dif" + accessionIndex + ".txt");
                        PrintWriter out3 = new PrintWriter(writer3);
                        for (int q=0; q<IndexesToGet.length; q++) {
                                out3.print(differences[q]);
                                out3.print("\t");
                        }
                        out3.print("\n");
                        out3.close();




	

			// Write accession names
			if (accessionIndex == 5) {
				Writer writerN = new FileWriter(results + "names.txt");
				PrintWriter outN = new PrintWriter(writerN);
	    	   		for (int sub=0; sub<IndexesToGet.length; sub++) {
					outN.print(splitIDs[IndexesToGet[sub]] + "\n");
				}
	       			outN.print("\n");
				outN.close();
			}
	
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		Pairwise_shore_clean pairwise_shore_clean = new Pairwise_shore_clean();
		pairwise_shore_clean.setMatrixFile(args[0], args[1], args[2], args[3]);
	}
}
