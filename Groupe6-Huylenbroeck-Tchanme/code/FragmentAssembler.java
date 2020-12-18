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
	*/ 
	private static int[][][] overlap_multigraph;


	/***/
	public static void main(String[] args){
		String path_in = "";
		String path_out = "";
		String path_ic_out = "";
		if(args.length==5){
			path_in = args[0];
			path_out = args[2];
			path_ic_out = args[4];
		}
		String collection_number = path_in.substring(path_in.length()-7, path_in.length()-6);
		String target_in = new File(path_in).getAbsoluteFile().getParent()+"/cible"+collection_number+".fasta";
		int target_len = OpenFasta(target_in).get(0).length;
		System.out.println("Fragment assembler by Huylenbroeck Florent and Tchanme Paul");
		fragments = OpenFasta(args[0]);
		CleanFragments();
		overlap_multigraph = OverlapMultigraph(fragments);
		int[] path = GreedyHamiltonianPath();
		byte[] contig = BuildConting(path);

		WriteFasta(path_out, contig, collection_number, Integer.toString(target_len));
		WriteFasta(path_ic_out, InvertAndComplement(contig), collection_number, Integer.toString(target_len));
		System.out.println("Done");
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
				if(line.length()>0){
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
			}
			if (fragment!="")
				fragments.add(fragment.getBytes());
		} catch(Exception e) {
			e.printStackTrace();
		}

		return fragments;
	}

	/**
	* Writes a .fasta file given an input sequence
	*
	* @param path 				String, path to the .fasta file.
	* @param sequence 			byte[], the sequence to be written.
	* @param collection_number 	String, the collection number from which the fragments have been read.
	* @param target_length 	 	String, length of the sequence from which the fragments are from.
	*/
	private static void WriteFasta(String path, byte[] sequence, String collection_number, String target_length){
		String fasta = "> Groupe-6 Collection "+collection_number+" "+sequence.length+" "+target_length+"\n";
		for(int i=0; i<sequence.length; i++){
			fasta = fasta.concat(Character.toString((char)sequence[i]));
		}
		try{
			FileWriter fastaWriter = new FileWriter(path);
			fastaWriter.write(fasta);
			fastaWriter.close();
		} catch(Exception e) {
			e.printStackTrace();
		}

		
	}


	/**
	* Removes fragments that are completely contained in another fragment, or it's inverted complement.
	*/
	private static void CleanFragments(){
		System.out.println("Removing contained fragments.");
		for(int i=0; i<fragments.size(); i++){
			for(int j=i; j<fragments.size(); j++){
				if(i!=j && (Contained(fragments.get(i), fragments.get(j)) || Contained(fragments.get(i), InvertAndComplement(fragments.get(j))))) {
					fragments.remove(fragments.get(j));
				}
			}
		}
	}

	/**
	* Tells if fragment 2 is contained in fragment 1
	*
	* @param fragment1 	byte[], the container fragment.
	* @param fragment2 	byte[], the potentially contained fragment.
	*/
	private static boolean Contained(byte[] fragment1, byte[] fragment2){
		int index1=0;
		int index2=0;
		while(index1<fragment1.length && index2<fragment2.length){
			if(fragment1[index1]==fragment2[index2]){
				index1++;
				index2++;
			}
			else{
				index1++;
			}
		}

		return (index2==fragment2.length);
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
				case (byte)97 : fragmentCI[fragment.length-i-1]=(byte)116; break;
				case (byte)99 : fragmentCI[fragment.length-i-1]=(byte)103; break;
				case (byte)103 : fragmentCI[fragment.length-i-1]=(byte)99; break;
				case (byte)116 : fragmentCI[fragment.length-i-1]=(byte)97; break;
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
	* @return 			int, semi-global alignment score.
	*/
	private static int SemiGlobalAlignmentScore(byte[] fragment1, byte[] fragment2){
		int m = fragment1.length, n = fragment2.length, gap_score = -2, mismatch_score = -1, match_score = 1, treshold = 3;

		int[] a = new int[n+1];

		for(int j = 0; j<=n; j++){
			a[j]= 0; //-1*(int)(j/treshold); //j*gap_score; // 0 to favor gaps in the beginning of the sequence
		}

		for(int i = 1; i<=m; i++){
			int old = a[0];
			a[0] = 0; //-1*(int)(i/treshold); // 0 to favor gaps in the beginning of the sequence
			for(int j = 1; j<=n; j++){
				int temp = a[j];
				int p = (fragment1[i-1] == fragment2[j-1]) ? match_score:mismatch_score; // -1 to shift the entries in the sims tab
				a[j] =  Math.max(a[j]+gap_score, Math.max(old+p, a[j-1]+gap_score));
				old = temp;
			}
		}

		return a[n];
	}

	/**
	* Computes the semi global alignment of two fragments.
	*
	* @param fragment1 	byte[] that holds the first fragment to be aligned.
	* @param fragment2 	byte[] that holds the second fragment to be aligned.
	* @return 			int[][] matrix that holds the alignment scores.
	*/
	private static int[][] BuildAlignmentMatrix(byte[] fragment1, byte[] fragment2){
		int m = fragment1.length+1, n = fragment2.length+1, gap_score = -2, mismatch_score = -1, match_score = 1, treshold = 3;;

		int[][] sims = new int[m][n];

		for(int i = 0; i<m; i++){
			sims[i][0] = 0;//-1*(int)(i/treshold);; // Favors gaps in the beginning of the sequence, was i*gap_score;
		}

		for(int j = 0; j<n; j++){
			sims[0][j] = 0;//-1*(int)(j/treshold);; // Favors gaps in the beginning of the sequence, was j*gap_score;
		}

		for(int i = 1; i<m; i++){
			for(int j = 1; j<n; j++){
				int p = (fragment1[i-1] == fragment2[j-1]) ? match_score:mismatch_score;
				sims[i][j] = Math.max(sims[i-1][j] + gap_score, Math.max(sims[i-1][j-1] + p, sims[i][j-1] + gap_score));
			}
		}

		return sims;
	}

	/**
	* Builds an overlap multigraph based on a set of fragments.
	*
	* @param fragments 	ArrayList<byte[]>, the fragments from which the graph will be built.
	* @return 			int[][][], the overlap multigraph.
	*/
	private static int[][][] OverlapMultigraph(ArrayList<byte[]> fragments){
		System.out.println("Building overlap multigraph.");
		overlap_multigraph = new int[fragments.size()][fragments.size()][4];

		for(int i=0; i<fragments.size(); i++){

			byte[] fragment1 = fragments.get(i);
			byte[] fragmentIC1 = InvertAndComplement(fragment1); //TODO reduce the number of call to InvertAndComplement()
			for(int j=0; j<fragments.size(); j++){
				if (i==j)
					continue;
				byte[] fragment2 = fragments.get(j);
				byte[] fragmentIC2 = InvertAndComplement(fragment2);

				overlap_multigraph[i][j][0] = SemiGlobalAlignmentScore(fragment1, fragment2); // 1 - 2 == 2(ic) - 1(ic)
				overlap_multigraph[i][j][1] = SemiGlobalAlignmentScore(fragment1, fragmentIC2); // 1 - 2(ic) == 2 - 1(ic)
				overlap_multigraph[i][j][2] = SemiGlobalAlignmentScore(fragmentIC1, fragment2); // 1(ic) - 2 == 2(ic) - 1
				overlap_multigraph[i][j][3] = SemiGlobalAlignmentScore(fragmentIC1, fragmentIC2); // 1(ic) - 2(ic) == 2 - 1
			}
		}
		return overlap_multigraph;
	}

	/**
	* Finds an hamiltonian path in an overlap multigraph
	*
	* @return 	int[] hamiltonian path. negative indexes correspond to inverted and complemented fragments.
	*/
	private static int[] GreedyHamiltonianPath(){
		System.out.println("Computing the hamiltonian path in the overlap multigraph.");

		// init
		byte[] in = new byte[overlap_multigraph.length];
		byte[] out = new byte[overlap_multigraph.length];
		ArrayList<int[]> sets = new ArrayList<int[]>();

		for(int i=0; i<overlap_multigraph.length; i++){
			sets.add(new int[] {i});
		}

		// vertice sorting in descending order.
		ArrayList<int[]> vertice = new ArrayList<int[]>();

		for(int i=0; i<overlap_multigraph.length; i++){
			for(int j=0; j<overlap_multigraph[i].length; j++){
				for(int k=0; k<4; k++){
					vertice.add(new int[] {overlap_multigraph[i][j][k], i, j, k});
				}
			}
		}

		vertice.sort(Comparator.comparing(a -> a[0]));
		Collections.reverse(vertice);
		//vertice.stream().map(Arrays::toString).forEach(System.out::println);

		int[][] greedy_hamiltionian_path_vertice = new int[overlap_multigraph.length-1][2];
		int greedy_index = 0;

		// pricess
		for(int i=0; i<vertice.size(); i++){
			int[] vertex = vertice.get(i);
			int[] set1 = FindSet(sets, vertex[1]);
			int[] set2 = FindSet(sets, vertex[2]);


			// test edge for acceptance
			if (in[vertex[1]] == 0 && out[vertex[2]] == 0 && !set1.equals(set2)){

				// ensure fragment are considered in their inverted and complemented state if necessary
				if (out[vertex[1]] > 0 && (vertex[3]==2 || vertex[3]==3)){
					continue;
				}
				if (out[vertex[1]] < 0 && (vertex[3]==0 || vertex[3]==1)){
					continue;
				}
				if (in[vertex[2]] > 0 && (vertex[3]==1 || vertex[3]==3)){
					continue;
				}
				if (in[vertex[2]] < 0 && (vertex[3]==0 || vertex[3]==2)){
					continue;
				}


				// ensure framgent will be considered in their inverted and complemented state if necessary
				if (vertex[3] == 0){
					in[vertex[1]] = 1;
					out[vertex[2]] = 1;
					greedy_hamiltionian_path_vertice[greedy_index]=new int[] {vertex[1], vertex[2]};
				}
				else if (vertex[3] == 1){
					in[vertex[1]] = 1;
					out[vertex[2]] = -1;
					greedy_hamiltionian_path_vertice[greedy_index]=new int[] {vertex[1], -vertex[2]};
				}
				else if (vertex[3] == 2){
					in[vertex[1]] = -1;
					out[vertex[2]] = 1;
					greedy_hamiltionian_path_vertice[greedy_index]=new int[] {-vertex[1], vertex[2]};
				}
				else{
					in[vertex[1]] = -1;
					out[vertex[2]] = -1;
					greedy_hamiltionian_path_vertice[greedy_index]=new int[] {-vertex[1], -vertex[2]};
				}
				greedy_index++;
				Union(sets, set1, set2);
			}
			if (sets.size()==1){
				break;
			}
		}

		// translate vertice list in a path
		int[] greedy_hamiltionian_path = new int[overlap_multigraph.length];

		// find starting point (the one that has no entry in 'in' but has one in 'out')
		for(int i=0; i<out.length; i++){
			if(in[i]==0){
				if(out[i]>0){
					greedy_hamiltionian_path[0]=i;
				}
				else{
					greedy_hamiltionian_path[0]=-i;
				}
				break;
			}
		}

		// then fill vector 
		for(int i=1; i<overlap_multigraph.length; i++){
			for(int[] vertex: greedy_hamiltionian_path_vertice){
				if(vertex[1]==greedy_hamiltionian_path[i-1]){
					greedy_hamiltionian_path[i]=vertex[0];
					break;
				}
			}
		}

		return greedy_hamiltionian_path;
	}

	private static int[] FindSet(ArrayList<int[]> sets, int f){
		for(int[] set:sets){
			for(int i=0; i<set.length; i++){
				if (set[i]==f){
					return set;
				}
			}
		}
		return new int[] {};
	}

	private static void Union(ArrayList<int[]> sets, int[] set1, int[] set2){
		int new_length = set1.length+set2.length;
		int[] new_set = new int[new_length];

		int fill_index=0;
		for(int i=0; i<set1.length; i++){
			new_set[fill_index]=set1[i];
			fill_index++;
		}
		for(int i=0; i<set2.length; i++){
			new_set[fill_index]=set2[i];
			fill_index++;
		}

		sets.remove(set1);
		sets.remove(set2);

		sets.add(new_set);
	}

	/**
	* Builds the contig according to the hamiltonian path found amongst the fragments
	* 
	* @param path 	int[] hamiltonian path in the overlap multigraph.
	* @return 		byte[] consensus sequence.
	*/
	private static byte[] BuildConting(int[] path){
		System.out.println("Building consensus contig");
		byte[] consensus = new byte[5000];
		byte[] fragment1 = fragments.get(Math.abs(path[1]));

		ArrayList<Byte>[] alignments = new ArrayList[fragments.size()];

		for(int f=0; f<path.length; f++){	
			byte[] fragment2 = fragments.get(Math.abs(path[f]));
			if(path[f]<0){
				fragment2 = InvertAndComplement(fragment2);
			}
			int[][] alignment_matrix = BuildAlignmentMatrix(fragment1, fragment2);
		
			int best_score=Integer.MIN_VALUE, read_x=0, read_y=0;

			for(int j=0; j<alignment_matrix.length; j++){
				if (alignment_matrix[j][alignment_matrix[j].length-1]>best_score){
					best_score = alignment_matrix[j][alignment_matrix[j].length-1];
					read_x=j;
					read_y=alignment_matrix[j].length-1;
				}
				else if(alignment_matrix[j][alignment_matrix[j].length-1]==best_score && alignment_matrix[j].length-1+j > read_x+read_y){
					best_score = alignment_matrix[j][alignment_matrix[j].length-1];
					read_x=j;
					read_y=alignment_matrix[j].length-1;
				}
				if(j==alignment_matrix.length-1){
					for(int k=0; k<alignment_matrix[j].length-1; k++){
						if(alignment_matrix[j][k]>best_score){
							best_score = alignment_matrix[j][k];
							read_x=j;
							read_y=k;
						}
						else if(alignment_matrix[j][k]==best_score && j+k > read_x+read_y){
							best_score = alignment_matrix[j][k];
							read_x=j;
							read_y=k;
						}
					}
				}
			}
			String tmp_contig="";

			ArrayList<Byte> aligned_fragment = new ArrayList<Byte>();

			if(read_x<alignment_matrix.length-1){
				for(int j=alignment_matrix.length-1; j>read_x; j--){
					tmp_contig=Character.toString((char)fragment1[j-1]).concat(tmp_contig);
					aligned_fragment.add(0, (byte)'-');
				}
			}

			if(read_y<alignment_matrix[alignment_matrix.length-1].length-1){
				for(int j=alignment_matrix[alignment_matrix.length-1].length-1; j>read_y; j--){
					tmp_contig=Character.toString((char)fragment2[j-1]).concat(tmp_contig);
					aligned_fragment.add(0, fragment2[j-1]);
					for(int k=0; k<f; k++){
						alignments[k].add((byte)'-');
					}
				}
			}

			while(read_x>0 && read_y>0){
				if(fragment1[read_x-1]==fragment2[read_y-1]){
					tmp_contig=Character.toString((char)fragment1[read_x-1]).concat(tmp_contig);
					aligned_fragment.add(0, fragment2[read_y-1]);
					read_y--;
					read_x--;
				}
				else{
					int upper = alignment_matrix[read_x-1][read_y];
					int upper_left = alignment_matrix[read_x-1][read_y-1];
					int right = alignment_matrix[read_x][read_y-1];

					int max = Math.max(Math.max(upper, upper_left), right);
					if(max == upper){
						tmp_contig="-".concat(tmp_contig);;
						aligned_fragment.add(0, (byte)'-');
						read_x--;
					}
					else if(max == right){
						tmp_contig="-".concat(tmp_contig);;
						aligned_fragment.add(0, fragment2[read_y-1]);
						for(int k=0; k<f; k++){
							alignments[k].add(read_x, (byte)'-');
						}
						read_y--;
					}
					else{
						tmp_contig="-".concat(tmp_contig);
						aligned_fragment.add(0, fragment2[read_y-1]);
						read_x--;
						read_y--;
					}
				}
			}
			while(read_x!=0){
				tmp_contig=Character.toString((char)fragment1[read_x-1]).concat(tmp_contig);
				aligned_fragment.add(0, (byte)'-');
				read_x--;
			}
			while(read_y!=0){
				tmp_contig=Character.toString((char)fragment2[read_y-1]).concat(tmp_contig);
				aligned_fragment.add(0, fragment2[read_y-1]);
				read_y--;
				for(int k=0; k<f; k++){
					alignments[k].add(0, (byte)'-');
				}
			}

			alignments[f]=aligned_fragment;

			fragment1=tmp_contig.getBytes();
			
			
			for(int k=0; k<fragment1.length; k++){
				if(fragment1[k]==(char)'-'){
					int count_A = 0;
					int count_C = 0;
					int count_G = 0;
					int count_T = 0;
					int count_Gaps = 0;
					for(int j=0; j<=f; j++){
						byte tmp = alignments[j].get(k);
						if(tmp==(byte)97)
							count_A++;
						else if(tmp==(byte)99)
							count_C++;
						else if(tmp==(byte)103)
							count_G++;
						else if(tmp==(byte)116)
							count_T++;
						else if(tmp==(byte)45)
							count_Gaps++;

					}
					int max = Math.max(Math.max(Math.max(count_A, count_C), count_T), count_G);

					if(count_Gaps==max && count_A+count_C+count_G+count_T==0)
						fragment1[k]=45;
					else if(count_A==max)
						fragment1[k]=97;
					else if(count_C==max)
						fragment1[k]=99;
					else if(count_G==max)
						fragment1[k]=103;
					else if(count_T==max)
						fragment1[k]=116;
				}
			}
			consensus=fragment1;
		}
		return consensus;
	}
}