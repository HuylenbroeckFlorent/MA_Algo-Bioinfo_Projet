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
	* Adjacency matrix for the overlap multigraph.
	* Dimensions are :
	* 1 : first fragment index.
	* 2 : second fragment index.
	* 3 : four possible edges between those two fragments, as follow :
	*		1 : Semiglobal alignment score between the both fragments.
	*		2 : Semiglobal alignment score between the first fragment and the inverted and complemented second fragment.
	*		3 : Semiglobal alignment score between the inverted and complemented first fragment and the second fragment.
	*		4 : Semiglobal alignment score between both inverted and complemented fragments.
	* This matrix is triangular, so that the same value is not stored twice for every edges.
	*/ 
	private static int[][][] overlap_multigraph;

	/***/
	public static void main(String[] args){
		fragments = OpenFasta(args[0]);
		overlap_multigraph = OverlapMultigraph(fragments);
		

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
	* (byte) '-' == 45
	* (byte) 'a' == 97
	* (byte) 'c' == 99
	* (byte) 'g' == 103
	* (byte) 't' == 116
	*
	* @param fragment 	byte[], fragment to be complemented and inverted.
	* @return 			byte[], the complemented and inverted fragment. 
	*/
	private static byte[] InvertAndComplement(byte[] fragment){
		byte[] fragmentCI = new byte[fragment.length];

		for(int i=0; i<fragment.length; i++){
			switch(fragment[i]){
				case 97 : fragmentCI[fragment.length-i-1]=116; break;
				case 99 : fragmentCI[fragment.length-i-1]=103; break;
				case 103 : fragmentCI[fragment.length-i-1]=99; break;
				case 116 : fragmentCI[fragment.length-i-1]=97; break;
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
	private static int SemiGlobalAlignmentScore(byte[] fragment1, byte[] fragment2){
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

		return sims[fragment1.length][fragment2.length];
	}

	/**
	* Builds an overlap multigraph based on a set of fragments.
	*
	* @param fragments 	ArrayList<byte[]>, the fragments from which the graph will be built.
	* @return 			int[][][], the overlap multigraph.
	*/
	private static int[][][] OverlapMultigraph(ArrayList<byte[]> fragments){
		overlap_multigraph = new int[fragments.size()][fragments.size()][4];

		for(int i=0; i<fragments.size(); i++){
			for(int j=0; j<fragments.size(); j++){
				if (i==j)
					continue;
				byte[] fragment1 = fragments.get(i);
				byte[] fragment2 = fragments.get(j);
				byte[] fragmentIC1 = InvertAndComplement(fragment1); //TODO reduce the number of call to InvertAndComplement()
				byte[] fragmentIC2 = InvertAndComplement(fragment2);

				overlap_multigraph[i][j][0] = SemiGlobalAlignmentScore(fragment1, fragment2); // 1 - 2
				overlap_multigraph[i][j][1] = SemiGlobalAlignmentScore(fragment1, fragmentIC2); // 1 - 2(ic)
				overlap_multigraph[i][j][2] = SemiGlobalAlignmentScore(fragmentIC1, fragment2); // 1(ic) - 2
				overlap_multigraph[i][j][3] = SemiGlobalAlignmentScore(fragmentIC1, fragmentIC2); // 1(ic) - 2(ic)
			}
		}
		return overlap_multigraph;
	}

	/**
	*/
	private static int[] GreedyHamiltonianPath(int[][][] overlap_multigraph){
		boolean[] in = new boolean[overlap_multigraph[0].length];
		boolean[] out = new boolean[overlap_multigraph[0].length];
		return new int[1]; 
 	}
}