package utility;

import org.kohsuke.args4j.Option;

public class CmdArgs
{

	@Option(name = "-model", usage = "Specify model")
	public String model = "KIDSDMM";

	@Option(name = "-dataname", usage = "Data name")
	public String dataname = "SearchSnippets";

	@Option(name = "-corpus", usage = "Specify path to topic modeling corpus")
	public String corpus = "dataset/SearchSnippets.txt";

	@Option(name = "-output", usage = "Specify path to save the results")
	public String output_dir = "results/";

	@Option(name = "-schema", usage = "Specify path to the file containing similarity")
	public String schema = "";

	@Option(name = "-threshold", usage = "Specify threshold")
	public double threshold = 0.5;

	@Option(name = "-weight", usage = "Specify weight in GPUDMM")
	public double weight = 0.3;

	@Option(name = "-filterSize", usage = "Filter less semantic related word pairs in GPUDMM")
	public int filterSize = 20;

	@Option(name = "-ntopics", usage = "Specify number of topics")
	public int ntopics = 300;

	@Option(name = "-alpha", usage = "Specify alpha")
	public double alpha = 0.1;

	@Option(name = "-beta", usage = "Specify beta")
	public double beta = 0.1;

	@Option(name = "-lambda", usage = "Specify mixture weight lambda")
	public int lambda = 3;

	@Option(name = "-niters", usage = "Specify number of iterations")
	public int niters = 1000;

	@Option(name = "-nBurnIn", usage = "Specify number of burnIn")
	public int nBurnIn = 500;

	@Option(name = "-twords", usage = "Specify number of top topical words")
	public int twords = 20;

	@Option(name = "-name", usage = "Specify a name to topic modeling experiment")
	public String expModelName = "";

	@Option(name = "-dir")
	public String dir = "";

	@Option(name = "-label")
	public String labelFile = "";

	@Option(name = "-prob")
	public String prob = "";

	@Option(name = "-paras", usage = "Specify path to hyper-parameter file")
	public String paras = "";

	@Option(name = "-topWordsDir")
	public String topWordsPath = "";

	@Option(name = "-vocabFile")
	public String vocabFile = "";

	@Option(name = "-wikiDir")
	public String wikiDir = "";

	@Option(name = "-topicCoherEval")
	public String topicCoherEval = "pmi";

	@Option(name = "-topTC")
	public int topTC = -1;

}
