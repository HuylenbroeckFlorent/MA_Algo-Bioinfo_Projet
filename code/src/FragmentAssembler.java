import java.io.*;
import java.util.*;

class FragmentAssembler{

	/**
	* ArrayList that holds the fragments.
	*/
	public static ArrayList<String> fragments;

	public static void main(String[] args){
		fragments = OpenFasta(args[0]);

	}


	/**
	* Opens a .fasta file and imports the fragments that it describes.
	*/
	public static ArrayList<String> OpenFasta(String path){
		ArrayList<String> fragments = new ArrayList<String>();
		try{
			BufferedReader fastareader = new BufferedReader(new FileReader(path));
			String line;
			String fragment = "";
			while((line = fastareader.readLine()) != null){
				if (line.charAt(0)=='>'){
					if (fragment!=""){
						fragments.add(fragment);
						fragment="";
					}
				}
				else{
					fragment+=line;
				}
			}
			if (fragment!="")
				fragments.add(fragment);
		} catch(Exception e) {
			e.printStackTrace();
		}

		return fragments;
	}

	/**
	* Finds every occurences of a word in a sequence (the naive way).
	*/
	public static int[] FindOccurencesNaive(String seq, String word){
		int comp = 0;
		if(seq.length()<word.length())
			return null;
		for(int i=0; i<seq.length()-(word.length()+1); i++){
			int j=0;
			while(j<word.length()-1){
				if(word.charAt(j)==seq.charAt(i+j)){
					j++;
					comp++;
				}
				else
					break;
			}
			if(j==word.length()-1)
				System.out.println(i);
		}
		System.out.println(comp);
		return null;
	}
}