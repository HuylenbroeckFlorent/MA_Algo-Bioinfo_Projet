import java.io.*;
import java.util.*;

class FragmentAssembler{

	/* NOTES
	* DOTMATCHER INVERSE LES AXES ???
	*
	* A><T
	* C><G
	*
	* (byte) '-' == 45
	* (byte) 'a' == 97
	* (byte) 'c' == 99
	* (byte) 'g' == 103
	* (byte) 't' == 116
	*/

	/**
	* ArrayList that holds the fragments.
	*/
	public static ArrayList<byte[]> fragments;

	/***/
	public static void main(String[] args){;
		fragments = OpenFasta(args[0]);

		SemiGlobalAlignment(new byte[]{(byte) 'a', (byte) 'a',(byte) 'a',(byte) 'c'},  new byte[]{(byte) 'a',(byte) 'g',(byte) 'c'});
	}


	/**
	* Opens a .fasta file and imports the fragments that it describes.
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
	* Computes the semi global alignment of two fragments.
	*
	*
	*
	*/
	private static int[][] SemiGlobalAlignment(byte[] fragment1, byte[] fragment2){
		int m = fragment1.length+1, n = fragment2.length+1, gap_score = -2, mismatch_score = -1, match_score = 1;

		int[][] sims = new int[m+1][n+1];

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

		for(int i = 1; i<m; i++){
			for(int j = 1; j<n; j++){
				System.out.print(sims[i][j]+" ");
			}
			System.out.println();
		}

		return sims;
	}
}