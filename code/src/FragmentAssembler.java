import java.io.*;
import java.util.*;

class FragmentAssembler{

	/* NOTES
	* DOTMATCHER INVERSE LES AXES ???
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
}