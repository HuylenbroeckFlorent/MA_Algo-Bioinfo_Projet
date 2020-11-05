import java.io.*;
import java.util.*;

class FragmentAssembler{

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
}