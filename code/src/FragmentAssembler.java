import java.io.*;
import java.util.*;

class FragmentAssembler{

	/* NOTES
	* DOTMATCHER INVERSE LES AXES ???
	*/

	/**
	* ArrayList that holds the fragments.
	*/
	private static ArrayList<byte[]> fragments;

	/**
	* Adjacency matrix for the overlap graph.
	*/
	private static int[][] overlap_graph;

	/***/
	public static void main(String[] args){
		fragments = OpenFasta(args[0]);
		overlap_graph = new int[fragments.size()][fragments.size()];

		for(int i=0; i<fragments.size(); i++){
			for(int j=0; j<fragments.size(); j++){
				if (i==j)
					continue;
				overlap_graph[i][j]=SemiGlobalAlignment(fragments.get(i), fragments.get(j))[fragments.get(i).length][fragments.get(j).length];
			}
		}
		for(int i=0; i<fragments.size(); i++){
			for(int j=0; j<i; j++){
				System.out.print(overlap_graph[i][j]+ " ");
			}
			System.out.println();
		}

	}


	/**
	* Opens a .fasta file and imports the fragments that it describes.
	* (byte) '-' == 45
	* (byte) 'a' == 97
	* (byte) 'c' == 99
	* (byte) 'g' == 103
	* (byte) 't' == 116
	*
	* @param path 	String, path to the .fasta file.
	* @return 		Arraylist<byte[]> that contains all the retreived fragments, stored as individuals byte arrays (one byte = one char) in an Arraylist.
	*/
	private static ArrayList<byte[]> OpenFasta(String path){
		ArrayList<byte[]> fragments = new ArrayList<byte[]>();
		try{
			BufferedReader fastaReader = new BufferedReader(new FileReader(path));
			String line;
			String fragment = "";
			while((line = fastaReader.readLine()) != null){
				if (line.charAt(0)=='>'){
					if (fragment!=""){
						fragments.add(fragment.getBytes());
						fragment="";
					}
				}
				else{
					fragment+=line;
				}
			}
			if (fragment!="")
				fragments.add(fragment.getBytes());
		} catch(Exception e) {
			e.printStackTrace();
		}

		return fragments;
	}


	/**
	* Complements and invert a fragment.
	* A <-> T
	* C <-> G
	*
	* @param fragment 	byte[], fragment to be complemented and inverted.
	* @return 			byte[], the complemented and inverted fragment. 
	*/
	private static byte[] InvertAndComplement(byte[] fragment){
		byte[] fragmentCI = new byte[fragment.length];

		for(int i=0; i<fragment.length; i++){
			switch(fragment[i]){
				case 'a' : fragmentCI[fragment.length-i-1]='t'; break;
				case 'c' : fragmentCI[fragment.length-i-1]='g'; break;
				case 'g' : fragmentCI[fragment.length-i-1]='c'; break;
				case 't' : fragmentCI[fragment.length-i-1]='a'; break;
				default : break;
			}
		}

		return fragmentCI;
	}


	/**
	* Computes the semi global alignment of two fragments.
	*
	* @param fragment1 	byte[] that holds the first fragment to be aligned.
	* @param fragment2 	byte[] that holds the second fragment to be aligned.
	* @return 			int[][] matrix that holds the alignment scores.
	*/
	private static int[][] SemiGlobalAlignment(byte[] fragment1, byte[] fragment2){
		int m = fragment1.length+1, n = fragment2.length+1, gap_score = -2, mismatch_score = -1, match_score = 1;

		int[][] sims = new int[m][n];

		for(int i = 0; i<m; i++){
			sims[i][0] = 0; // Favors gaps in the beginning of the sequence, was i*gap_score;
		}

		for(int j = 0; j<n; j++){
			sims[0][j] = 0; // Favors gaps in the beginning of the sequence, was j*gap_score;
		}

		for(int i = 1; i<m; i++){
			for(int j = 1; j<n; j++){
				int p = (fragment1[i-1] == fragment2[j-1]) ? 1:-1; // -1 to shift the entries in the sims tab
				sims[i][j] = Math.max(sims[i-1][j] + gap_score, Math.max(sims[i-1][j-1] + p, sims[i][j-1] + gap_score));
			}
		}

		return sims;
	}
}