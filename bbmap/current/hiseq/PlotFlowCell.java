package hiseq;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import jgi.Dedupe;
import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.ScheduleMaker;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Analyzes a flow cell for low-quality areas.
 * Removes reads in the low-quality areas.
 * 
 * @author Brian Bushnell
 * @date August 31, 2016
 *
 */
public class PlotFlowCell {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		PlotFlowCell x=new PlotFlowCell(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public PlotFlowCell(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		//Create a parser object
		Parser parser=new Parser();
		parser.qtrimRight=trimRight;
		parser.trimq=trimq;
		parser.minReadLength=minlen;
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("divisor") || a.equals("size")){
				Tile.xSize=Tile.ySize=Integer.parseInt(b);
			}else if(a.equals("xdivisor") || a.equals("xsize")){
				Tile.xSize=Integer.parseInt(b);
			}else if(a.equals("ydivisor") || a.equals("ysize")){
				Tile.ySize=Integer.parseInt(b);
			}else if(a.equals("target")){
				targetAverageReads=Integer.parseInt(b);
			}else if(a.equals("dump")){
				dump=b;
			}else if(a.equals("indump") || a.equals("ind") || a.equals("dumpin")){
				dumpIn=b;
			}else if(a.equals("pound")){
				pound=Tools.parseBoolean(b);
			}else if(a.equals("loadkmers") || a.equals("usekmers")){
				loadKmers=Tools.parseBoolean(b);
			}else if(a.equals("lqo") || a.equals("lowqualityonly")){
				discardOnlyLowQuality=Tools.parseBoolean(b);
			}else if(a.equals("dl") || a.equals("discardlevel")){
				discardLevel=Integer.parseInt(b);
			}else if(a.equals("outbad") || a.equals("outb") || a.equals("outtoss") || a.equals("outt") || a.equals("outunwanted") || a.equals("outu")){
				outbad=b;
			}
			
			else if(a.equals("deviations") || a.equals("d")){
				qDeviations=uDeviations=eDeviations=Float.parseFloat(b);
			}else if(a.equals("qdeviations") || a.equals("qd") || a.equals("dq")){
				qDeviations=Float.parseFloat(b);
			}else if(a.equals("udeviations") || a.equals("ud") || a.equals("du")){
				uDeviations=Float.parseFloat(b);
			}else if(a.equals("edeviations") || a.equals("ed") || a.equals("de")){
				eDeviations=Float.parseFloat(b);
			}else if(a.equals("qfraction") || a.equals("qf")){
				qualFraction=Float.parseFloat(b);
			}else if(a.equals("efraction") || a.equals("uf")){
				uniqueFraction=Float.parseFloat(b);
			}else if(a.equals("efraction") || a.equals("ef")){
				errorFreeFraction=Float.parseFloat(b);
			}else if(a.equals("qabsolute") || a.equals("qa")){
				qualAbs=Float.parseFloat(b);
			}else if(a.equals("uabsolute") || a.equals("ua")){
				uniqueAbs=Float.parseFloat(b);
			}else if(a.equals("eabsolute") || a.equals("ea")){
				errorFreeAbs=Float.parseFloat(b);
			}
			
			else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
			

			trimq=parser.trimq;
			trimE=parser.trimE();
			minlen=parser.minReadLength;
			trimLeft=parser.qtrimLeft;
			trimRight=parser.qtrimRight;
		}
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outbad, dump)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2+", "+outbad);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+", "+outbad+", "+dump+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outbad, dump)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public void process(Timer t){
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		if(dumpIn==null){
			flowcell=new FlowCell();
			if(loadKmers){loadKmers();}
			fillTiles();
			keySets=null;
		}else{
			flowcell=new FlowCell(dumpIn);
			
			if(flowcell.avgReads<targetAverageReads){
				flowcell=flowcell.widen(targetAverageReads);
			}
			
			avgQuality=flowcell.avgQuality;
			avgUnique=flowcell.avgUnique;
			avgErrorFree=flowcell.avgErrorFree;
			stdQuality=flowcell.stdQuality;
			stdUnique=flowcell.stdUnique;
			stdErrorFree=flowcell.stdErrorFree;
		}
	}

	/** Create read streams and process all data */
	void loadKmers(){
		Timer t2=new Timer();
		outstream.print("Loading kmers:  \t");
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Process the read stream
		loadKmersInner(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		t2.stop();
		outstream.println(t2);
	}

	/** Create read streams and process all data */
	void fillTiles(){
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Process the read stream
		fillTilesInner(cris);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
	}
	
	/** Iterate through the reads */
	public void loadKmersInner(final ConcurrentReadInputStream cris){
		
		keySets=new AbstractKmerTable[WAYS];

		//Initialize tables
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, false, 0.8);
		int[] schedule=scheduleMaker.makeSchedule();
		for(int j=0; j<WAYS; j++){
			keySets[j]=new HashArray1D(schedule, -1L);
		}
//		for(int j=0; j<WAYS; j++){
//			keySets[j]=new HashArray1D(512000, -1, -1L, true); //TODO: Set maxSize
//		}
		//Do anything necessary prior to processing
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					if(initialLength1>=k && randy.nextBoolean()){
						final long kmer=toKmer(r1.bases, randy.nextInt(initialLength1-k2), k);
						if(kmer>=0){
							AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
							table.increment(kmer, 1);
						}
					}
					
					if(initialLength2>=k && randy.nextBoolean()){
						final long kmer=toKmer(r2.bases, randy.nextInt(initialLength2-k2), k);
						if(kmer>=0){
							AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
							table.increment(kmer, 1);
						}
					}
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		//Do anything necessary after processing
	}
	
	/** Iterate through the reads */
	public void fillTilesInner(final ConcurrentReadInputStream cris){


		Timer t2=new Timer();
		outstream.print("Filling tiles:  \t");
		
		
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					//Increment counters
					readsProcessed+=r1.pairCount();
					basesProcessed+=initialLength1+initialLength2;
					
					final MicroTile mt=flowcell.getMicroTile(r1.id);
					
					if(loadKmers){
						if(initialLength1>=k){
							final long kmer=toKmer(r1.bases, randy.nextInt(initialLength1-k2), k);
							if(kmer>=0){
								AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
								if(table.getValue(kmer)>0){mt.hits++;}
								else{mt.misses++;}
							}else{mt.misses++;}
						}
						
						if(initialLength1>=k){
							final long kmer=toKmer(r1.bases, randy.nextInt(initialLength1-k2), k);
							if(kmer>=0){
								AbstractKmerTable table=keySets[(int)(kmer%WAYS)];
								if(table.getValue(kmer)>0){mt.hits++;}
								else{mt.misses++;}
							}else{mt.misses++;}
						}
					}

					mt.add(r1);
					mt.add(r2);
				}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				
				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		t2.stop();
		outstream.println(t2);
		
		ArrayList<MicroTile> mtList=flowcell.calcStats();
		if(flowcell.avgReads<targetAverageReads){
			flowcell=flowcell.widen(targetAverageReads);
			mtList=flowcell.toList();
		}
		avgQuality=flowcell.avgQuality;
		avgUnique=flowcell.avgUnique;
		avgErrorFree=flowcell.avgErrorFree;
		
		stdQuality=flowcell.stdQuality;
		stdUnique=flowcell.stdUnique;
		stdErrorFree=flowcell.stdErrorFree;
		
		if(dump!=null){
			TextStreamWriter tsw=new TextStreamWriter(dump, overwrite, append, false);
			tsw.start();

			tsw.println("#xSize\t"+Tile.xSize);
			tsw.println("#ySize\t"+Tile.ySize);
			tsw.println("#reads\t"+String.format(Locale.ROOT, "%d", flowcell.readsProcessed));
			tsw.println("#avgReads\t"+String.format(Locale.ROOT, "%.1f", flowcell.avgReads));
			
			tsw.println("#avgQuality\t"+String.format(Locale.ROOT, "%.3f", avgQuality));
			tsw.println("#avgUnique\t"+String.format(Locale.ROOT, "%.3f", avgUnique));
			tsw.println("#avgErrorFree\t"+String.format(Locale.ROOT, "%.3f", avgErrorFree));
			
			tsw.println("#stdQuality\t"+String.format(Locale.ROOT, "%.5f", stdQuality));
			tsw.println("#stdUnique\t"+String.format(Locale.ROOT, "%.5f", stdUnique));
			tsw.println("#stdErrorFree\t"+String.format(Locale.ROOT, "%.5f", stdErrorFree));
			
			tsw.println((pound ? "#" : "") + "lane\ttile\tx1\tx2\ty1\ty2\treads\tunique\tquality\tprobErrorFree\tdiscard\tA\tC\tG\tT\tN\tpolyG\tpolyGLen");
			
			for(Lane lane : flowcell.lanes){
				if(lane!=null){
					for(Tile tile : lane.tiles){
						if(tile!=null){
							tsw.print(tile.toString());
						}
					}
				}
			}
			tsw.poisonAndWait();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Generate a kmer from specified start location
	 * @param bases
	 * @param start
	 * @param klen kmer length
	 * @return kmer
	 */
	private static final long toKmer(final byte[] bases, final int start, final int klen){
		final int stop=start+klen;
		assert(stop<=bases.length) : klen+", "+bases.length;
		long kmer=0;
		
		for(int i=start; i<stop; i++){
			final byte b=bases[i];
			final long x=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|x);
		}
		return kmer;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	/** Discard output file path */
	private String outbad=null;

	private String qfout1=null;
	private String qfout2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	private boolean pound=true;
	private String dump=null;
	private String dumpIn=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	public long readsProcessed=0;
	/** Number of bases processed */
	public long basesProcessed=0;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private AbstractKmerTable[] keySets;
	
	private int targetAverageReads=800;

	private float minprob=0;
	private static final int WAYS=31;
	private static final int k=31, k2=30;
	
	private final Random randy=Shared.threadLocalRandom();
	private FlowCell flowcell;
	
//	private int[] avgQualityArray;
//	private int[] avgUniqueArray;
	
	private long minCountToUse=0;

	private float qDeviations=2f;
	private float uDeviations=1.5f;
	private float eDeviations=2f;
	
	private float qualFraction=0.01f;
	private float uniqueFraction=0.01f;
	private float errorFreeFraction=0.01f;
	
	private float qualAbs=1f;
	private float uniqueAbs=1f;
	private float errorFreeAbs=1f;
	
	private double avgQuality;
	private double avgUnique;
	private double avgErrorFree;
	
	private double stdQuality;
	private double stdUnique;
	private double stdErrorFree;

	private boolean loadKmers=true;
	private boolean discardOnlyLowQuality=true;
	private int discardLevel=1;
	
	private int minlen=30;
	private float trimq=-1;
	private final float trimE;
	private boolean trimLeft=false;
	private boolean trimRight=true;
	
	private boolean warned=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	/** Secondary input file */
	private final FileFormat ffin2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}
