package models;

import utility.FuncUtils;

import java.io.*;
import java.util.*;

public class KIDSDMM {

    public double alpha, beta;
    public int lambda;
    public int numTopics; // Number of topics
    public int numIterations; // Number of Gibbs sampling iterations
    public int numBurnIn; // Number of BurnIn
    public int top; // Number of most probable words for each topic
    public double alphaSum; // alpha * numTopics
    public double betaSum; // beta * vocabularySize

    public ArrayList<int[]> Corpus = new ArrayList<>(); // Word ID-based corpus

    private Random rg;
    public double threshold;

    public int filterSize;

    public HashMap<String, Integer> word2IdVocabulary; // Vocabulary to get ID
    // given a word
    public HashMap<Integer, String> id2WordVocabulary; // Vocabulary to get word
    // given an ID
    // Given a document, number of times its i^{th} word appearing from
    // the first index to the i^{th}-index in the document
    // Example: given a document of "a a b a b c d c". We have: 1 2 1 3 2 1 1 2
    public List<List<Integer>> occurenceToIndexCount;

    public int vocabularySize; // The number of word types in the corpus

    public int numDocuments; // Number of documents in the corpus
    public int numWordsInCorpus; // Number of words in the corpus

    public double[][] phi;
    private double[] pz;
    private double[][] pdz;
    private double[][] topicProbabilityGivenWord;

    public ArrayList<ArrayList<Boolean>> wordGPUFlag; // wordGPUFlag.get(doc).get(wordIndex)
    public int[] assignmentList; // topic assignment for every document
    public ArrayList<ArrayList<Map<Integer, Double>>> wordGPUInfo;

    // Number of documents assigned to a topic
    public int[] docTopicCount;
    // numTopics * vocabularySize matrix
    // Given a topic: number of times a word type assigned to the topic
    public double[][] topicWordCount;
    // Total number of words assigned to a topic
    public double[] sumTopicWordCount;

    private Map<Integer, Map<Integer, Double>> schemaMap;
    private double[][] schema;
    //  Local word correlation
    private double[][] LC;

    // Double array used to sample a topic
    public double[] multiPros;

    // Path to the directory containing the corpus
    public String folderPath;
    // Path to the topic modeling corpus
    public String corpusPath;

    public String schemafile;
    public double count;

    public String expName = "KIDSDMM model";
    public String orgExpName = "KIDSDMM model";
    public String tAssignsFilePath = "";

    public double[][] TWC; // V * K, TWC[v][k] : the correlation between word v and topic k
    public int[][] topWords; // K * top, topWords[k] : the index of most probable words in topic k

    public double[] sigma;

    public double[] oldWordWeight;
    public double[] bdcWeight; // bdcWeight[i] : the weight for word i

    public int savestep = 0;

    public double initTime = 0;
    public double iterTime = 0;

    public KIDSDMM(String pathToCorpus, String pathToResult, String schemaPath, double inThreshold,
                   int inFilterSize, int inNumTopics, double inAlpha, double inBeta, int inLambda, int inNumIterations,
                   int inNumBurnIn, int inTopWords, String inExpName)
            throws Exception {
        this(pathToCorpus, pathToResult, schemaPath, inThreshold,
                inFilterSize, inNumTopics, inAlpha, inBeta, inLambda, inNumIterations,
                inNumBurnIn, inTopWords, inExpName, 0);
    }

    public KIDSDMM(String pathToCorpus, String pathToResult, String schemaPath, double inThreshold,
                   int inFilterSize, int inNumTopics, double inAlpha, double inBeta, int inLambda, int inNumIterations,
                   int inNumBurnIn, int inTopWords, String inExpName, int inSaveStep)
            throws Exception {

        alpha = inAlpha;
        beta = inBeta;
        lambda = inLambda;
        numTopics = inNumTopics;
        numIterations = inNumIterations;
        numBurnIn = inNumBurnIn;
        top = inTopWords;
        expName = inExpName;
        orgExpName = expName;
        filterSize = inFilterSize;
        threshold = inThreshold;
        corpusPath = pathToCorpus;
        savestep = inSaveStep;

        schemafile = schemaPath;
        count = 0.0;
        folderPath = pathToResult;
        File dir = new File(folderPath);
        if (!dir.exists())
            dir.mkdirs();

        word2IdVocabulary = new HashMap<String, Integer>();
        id2WordVocabulary = new HashMap<Integer, String>();
        occurenceToIndexCount = new ArrayList<List<Integer>>();

        wordGPUFlag = new ArrayList<>();
        numDocuments = 0;
        numWordsInCorpus = 0;

        rg = new Random();

        BufferedReader br = null;
        try {
            int indexWord = -1;
            br = new BufferedReader(new FileReader(pathToCorpus));
            for (String doc; (doc = br.readLine()) != null;) {

                if (doc.trim().length() == 0)
                    continue;
                String[] words = doc.trim().split("\\s+");
                int [] document = new int[words.length];

                List<Integer> wordOccurenceToIndexInDoc = new ArrayList<Integer>();
                HashMap<Integer, Integer> wordOccurenceToIndexInDocCount = new HashMap<Integer, Integer>();

                int ind = 0;
                for (String word : words) {
                    if (word2IdVocabulary.containsKey(word)) {
                        document[ind++] = word2IdVocabulary.get(word);
                    }
                    else {
                        indexWord += 1;
                        word2IdVocabulary.put(word, indexWord);
                        id2WordVocabulary.put(indexWord, word);
                        document[ind++] = indexWord;
                    }
                    int times = 0;
                    if (wordOccurenceToIndexInDocCount.containsKey(indexWord)) {
                        times = wordOccurenceToIndexInDocCount.get(indexWord);
                    }
                    times += 1;
                    wordOccurenceToIndexInDocCount.put(indexWord, times);
                    wordOccurenceToIndexInDoc.add(times);
                }

                numDocuments++;
                numWordsInCorpus += document.length;
                Corpus.add(document);
                occurenceToIndexCount.add(wordOccurenceToIndexInDoc);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        vocabularySize = word2IdVocabulary.size(); // vocabularySize = indexWord
        docTopicCount = new int[numTopics];
        topicWordCount = new double[numTopics][vocabularySize];
        sumTopicWordCount = new double[numTopics];

        phi = new double[numTopics][vocabularySize];
        pz = new double[numTopics];

        topicProbabilityGivenWord = new double[vocabularySize][numTopics];

        pdz = new double[numDocuments][numTopics];
        multiPros = new double[numTopics];
        for (int i = 0; i < numTopics; i++) {
            multiPros[i] = 1.0 / numTopics;
        }

        alphaSum = numTopics * alpha;
        betaSum = vocabularySize * beta;

        assignmentList = new int[numDocuments];
        wordGPUInfo = new ArrayList<>();
        rg = new Random();

        LC = new double[vocabularySize][vocabularySize];
        schema = new double[vocabularySize][vocabularySize];

        TWC = new double[vocabularySize][numTopics];
        topWords = new int[numTopics][top];
        sigma = new double[vocabularySize];

        oldWordWeight = new double[vocabularySize];
        bdcWeight = new double[vocabularySize];
        for (int i = 0; i < vocabularySize; i++) {
            sigma[i] = 1.;
            bdcWeight[i] = 1.;
            oldWordWeight[i] = bdcWeight[i];
        }

        long startTime = System.currentTimeMillis();

        schemaMap = loadSchema(schemafile, threshold);
        // Compute the local word correlation
        computeLC();

        initialize();
        initTime =System.currentTimeMillis()-startTime;
        System.out.println("Corpus size: " + numDocuments + " docs, "
                + numWordsInCorpus + " words");
        System.out.println("Vocabulary size: " + vocabularySize);
        System.out.println("Number of topics: " + numTopics);
        System.out.println("alpha: " + alpha);
        System.out.println("beta: " + beta);
        System.out.println("threshold: " + threshold);
        System.out.println("filterSize: " + filterSize);
        System.out.println("Number of sampling iterations: " + numIterations);
        System.out.println("Number of top topical words: " + top);
        System.out.println(expName);

    }

    /**
     * Compute the Local word correlation
     */
    private void computeLC() {
        // wordTF[i] : word i term frequency
        double[] wordTF = new double[vocabularySize];
        // wordCount[i][j] : the probability of co-occurrence of word i and word j
        double[][] wordCount = new double[vocabularySize][vocabularySize];

        for (int docIndex=0; docIndex < Corpus.size(); docIndex++){
            int[] doc = Corpus.get(docIndex);
            for (int i = 0; i < doc.length; i++) {
                wordTF[doc[i]] += 1.0;
                for (int j = 0; j < doc.length; j++) {
                    if (i != j){
                        wordCount[doc[i]][doc[j]] += 1.0;
                    }
                }
            }
        }
        double numWords = 0.0, numCorrelation = 0.0;
        for (int i = 0; i < vocabularySize; i++) {
            numWords += wordTF[i];
            for (int j = 0; j < vocabularySize; j++) {
                numCorrelation += wordCount[i][j];
            }
        }
        for (int i = 0; i < vocabularySize; i++) {
            wordTF[i] /= numWords;
            for (int j = 0; j < vocabularySize; j++) {
                wordCount[i][j] /= numCorrelation;
            }
        }
        // Compute local PMI
        double result = 0.0;
        for (int i = 0; i < vocabularySize; i++) {
            for (int j = 0; j < vocabularySize; j++) {
                if (wordCount[i][j] != 0.0){
                    result = Math.log(wordCount[i][j] / (wordTF[i] * wordTF[j]));
                    if (result < 0) {
                        result = 0.0;
                    }
                    LC[i][j] = result;
                }
            }
        }
        System.out.println("Finish to compute local word correlation!\n");
    }

    private void updateWordWeight(){
        for (int i = 0; i < vocabularySize; i++) {
            oldWordWeight[i] = bdcWeight[i];
        }
    }

    /**
     * update sigma(w)
     *      sigma(w) = exp{(lambda - OF(w)) / K}
     * @return
     */
    private void updateSigma() {
        // Count the num of word occurrence
        Map<Integer, Set<Integer>> cnt = new HashMap<>();
        for (int i = 0; i < topWords.length; i++) {
            for (int j = 0; j < topWords[0].length; j++) {
                int wordId = topWords[i][j];
                Set<Integer> topicSet = new HashSet<>();
                if (cnt.containsKey(wordId)){
                    topicSet = cnt.get(wordId);
                }
                topicSet.add(i);
                cnt.put(wordId, topicSet);
            }
        }
        Set<Integer> wordIds = cnt.keySet();
        // Compute the mean of TU of top word
        for (int wordId: wordIds) {
            Set<Integer> topicSet = cnt.get(wordId);
            sigma[wordId] = Math.exp((lambda - topicSet.size()) / numTopics);
        }
    }

    private void updateBDC() {
        for (int i = 0; i < vocabularySize; i++) {
            double sum = 0.0;
            for (int j = 0; j < numTopics; j++) {
                double pro = computeBDCPro(i, j);
                sum += pro * Math.log(pro);
                if (Double.isNaN(sum)) {
                    System.out.println("***************");
                }
            }
            bdcWeight[i] = (1 + sum / Math.log(numTopics)) * sigma[i];
        }
    }

    private double computeBDCPro(int wordId, int topicId) {
        double sum = 0;
        for (int i = 0; i < numTopics; i++) {
            sum += phi[i][wordId];
        }
        return phi[topicId][wordId] / sum;
    }

    /**
     * Collect the similar words Map, not including the word itself
     *
     * @param filename:
     *            shcema_similarity filename
     * @param threshold:
     *            if the similarity is bigger than threshold, we consider it as
     *            similar words
     * @return
     */
    public Map<Integer, Map<Integer, Double>> loadSchema(String filename, double threshold) {
        int word_size = word2IdVocabulary.size();
        Map<Integer, Map<Integer, Double>> schemaMap = new HashMap<Integer, Map<Integer, Double>>();
        try {
            System.out.println(filename);
            FileInputStream fis = new FileInputStream(filename);
            InputStreamReader isr = new InputStreamReader(fis, "UTF-8");
            BufferedReader reader = new BufferedReader(isr);
            String line;
            int lineIndex = 0;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                String[] items = line.split(" ");

                for (int i = 0; i < items.length; i++) {
                    double value = Double.parseDouble(items[i]);
                    schema[lineIndex][i] = value;
                }
                lineIndex++;
            }

            for (int i = 0; i < word_size; i++) {
                Map<Integer, Double> tmpMap = new HashMap<Integer, Double>();
                for (int j = 0; j < word_size; j++) {
                    double v = schema[j][i];
                    if (Double.compare(v, threshold) > 0) {
                        tmpMap.put(j, v);
                    }
                }
                if (tmpMap.size() > filterSize) {
                    tmpMap.clear();
                }
                tmpMap.remove(i);
                if (tmpMap.size() == 0) {
                    continue;
                }
                count += tmpMap.size();
                schemaMap.put(i, tmpMap);
            }
            System.out.println("schema size is " + schemaMap.size());
            System.out.println("finish read schema, the average number of value is " + count / schemaMap.size());
            return schemaMap;
        } catch (Exception e) {
            System.out.println("Error while reading other file:" + e.getMessage());
            e.printStackTrace();
            return null;
        }
    }

    /**
     * Randomly initialize topic assignments
     */
    public void initialize()
            throws IOException
    {
        System.out.println("Randomly initializing topic assignments ...");

        for (int d = 0; d < numDocuments; d++) {
            int termIDArray[] = Corpus.get(d);

            ArrayList<Boolean> docWordGPUFlag = new ArrayList<>();
            ArrayList<Map<Integer, Double>> docWordGPUInfo = new ArrayList<>();

            ArrayList<int[]> d_assignment_list = new ArrayList<int[]>();

            int topic = FuncUtils.nextDiscrete(multiPros); // Sample a topic
            assignmentList[d] = topic;
            docTopicCount[topic] += 1;
            for (int t = 0; t < termIDArray.length; t++) {
                int termID = termIDArray[t];

                topicWordCount[topic][termID] += bdcWeight[termID];
                sumTopicWordCount[topic] += bdcWeight[termID];

                docWordGPUFlag.add(false); // initial for False for every word
                docWordGPUInfo.add(new HashMap<Integer, Double>());

            }
            wordGPUFlag.add(docWordGPUFlag);
            wordGPUInfo.add(docWordGPUInfo);

        }
        System.out.println("finish init_GPU!");

    }

    /**
     * Update S, the similarity matrix between topic k and word v
     *
     */
    public void updateWordTopicSimilarity() {
        for(int w = 0 ; w < vocabularySize ; w++) {
            for(int z = 0;z < numTopics; z++) {
                double tSim = 0;
                for(int v : topWords[z]) {
                    tSim += phi[z][v] * (schema[w][v] + LC[w][v]);
                }
                TWC[w][z] = Math.exp(tSim * topicProbabilityGivenWord[w][z]);
            }
        }
    }

    public void compute_phi() {
        for (int i = 0; i < numTopics; i++) {
            for (int j = 0; j < vocabularySize; j++) {
                phi[i][j] = (topicWordCount[i][j] + beta) / (sumTopicWordCount[i] + beta*vocabularySize);
                if (Double.isNaN(phi[i][j])) {
                    System.out.println("*********");
                }
            }
        }
    }

    public void compute_pz() {
        double sum = 0.0;
        for (int i = 0; i < numTopics; i++) {
            sum += sumTopicWordCount[i];
        }
        for (int i = 0; i < numTopics; i++) {
            pz[i] = (sumTopicWordCount[i] + alpha) / (sum + alphaSum);
        }
    }

    public void compute_pzd() {
        double[][] pwz = new double[vocabularySize][numTopics]; // pwz[word][topic]
        for (int i = 0; i < vocabularySize; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < numTopics; j++) {
                pwz[i][j] = pz[j] * phi[j][i];
                row_sum += pwz[i][j];
            }
            for (int j = 0; j < numTopics; j++) {
                pwz[i][j] = pwz[i][j] / row_sum;
            }
        }

        for (int i = 0; i < numDocuments; i++) {
            int[] doc_word_id = Corpus.get(i);
            double row_sum = 0.0;
            for (int j = 0; j < numTopics; j++) {
                for (int wordID : doc_word_id) {
                    pdz[i][j] += pwz[wordID][j];
                }
                row_sum += pdz[i][j];

            }
            for (int j = 0; j < numTopics; j++) {
                pdz[i][j] = pdz[i][j] / row_sum;
            }
        }
    }
    /**
     * update the p(z|w) for every iteration
     */
    public void updateTopicProbabilityGivenWord() {
        // TODO we should update pz and phi information before
        compute_pz();
        compute_phi();  //update p(w|z)
        for (int i = 0; i < vocabularySize; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < numTopics; j++) {
                topicProbabilityGivenWord[i][j] = pz[j] * phi[j][i];
                row_sum += topicProbabilityGivenWord[i][j];
            }
            for (int j = 0; j < numTopics; j++) {
                topicProbabilityGivenWord[i][j] = topicProbabilityGivenWord[i][j] / row_sum;  //This is p(z|w)
            }
        }
    }

    /**
     * Update topic top topical words
     */
    public void updateTopicTopWord() {
        for (int i = 0; i < phi.length; i++) {
            double[] phi_t = phi[i].clone();
            // get the topWords in the i_th topic
            for (int t = 0; t < topWords[0].length; t++) {
                int maxIndex = getMaxIndex(phi_t);
                // update top topical words in the i_th topic
                topWords[i][t] = maxIndex;
                phi_t[maxIndex] = 0;
            }
        }
    }

    /**
     * Get the max index in the double array
     * @param array
     * @return
     */
    private int getMaxIndex(double[] array) {
        int maxIndex = 0;
        double max = array[0];
        for (int i = 1; i < array.length; i++) {
            if (array[i] > max){
                max = array[i];
                maxIndex = i;
            }
        }
        return maxIndex;
    }

    public double findTopicMaxProbabilityGivenWord(int wordID) {
        double max = -1.0;
        for (int i = 0; i < numTopics; i++) {
            double tmp = topicProbabilityGivenWord[wordID][i];
            if (Double.compare(tmp, max) > 0) {
                max = tmp;
            }
        }
        return max;
    }

    /**
     * update GPU flag, which decide whether do GPU operation or not
     * @param docID
     * @param newTopic
     */
    public void updateWordGPUFlag(int docID, int newTopic) {
        // we calculate the p(t|w) and p_max(t|w) and use the ratio to decide we
        // use gpu for the word under this topic or not
        int[] termIDArray = Corpus.get(docID);
        ArrayList<Boolean> docWordGPUFlag = new ArrayList<>();
        for (int t = 0; t < termIDArray.length; t++) {

            int termID = termIDArray[t];
            double maxProbability = findTopicMaxProbabilityGivenWord(termID);
            if (maxProbability > 0) {
                double ratio = topicProbabilityGivenWord[termID][newTopic] / maxProbability;
                double a = rg.nextDouble();
                docWordGPUFlag.add(Double.compare(ratio, a) > 0);
            } else{
                docWordGPUFlag.add(false);
            }

        }
        wordGPUFlag.set(docID, docWordGPUFlag);
    }

    public void ratioCounter(Integer topic, Integer docID, int[] termIDArray, int flag) {
        docTopicCount[topic] += flag;
        // we update gpu flag for every document before it change the counter
        // when adding numbers
        if (flag > 0) {
            for (int t = 0; t < termIDArray.length; t++) {
                int wordID = termIDArray[t];
                topicWordCount[topic][wordID] += flag * bdcWeight[wordID];
                sumTopicWordCount[topic] += flag * bdcWeight[wordID];
            }
            updateWordGPUFlag(docID, topic);
            for (int t = 0; t < termIDArray.length; t++) {
                int wordID = termIDArray[t];
                boolean gpuFlag = wordGPUFlag.get(docID).get(t);
                Map<Integer, Double> gpuInfo = new HashMap<>();
                if (gpuFlag) { // do gpu count
                    if (schemaMap.containsKey(wordID)) {
                        Map<Integer, Double> valueMap = schemaMap.get(wordID);
                        // update the counter
                        for (Map.Entry<Integer, Double> entry : valueMap.entrySet()) {
                            Integer wordIndex = entry.getKey();
                            // Use the similarity of topic and word to update
                            double v = TWC[wordIndex][topic];
                            topicWordCount[topic][wordIndex] += v * bdcWeight[wordIndex];
                            sumTopicWordCount[topic] += v * bdcWeight[wordIndex];
                            gpuInfo.put(wordIndex, v);
                        } // end loop for similar words
                    } else { // schemaMap don't contain the word
                        // the word doesn't have similar words, the infoMap is empty
                        // we do nothing
                    }
                } else { // the gpuFlag is False
                    // it means we don't do gpu, so the gouInfo map is empty
                }
                wordGPUInfo.get(docID).set(t, gpuInfo); // we update the gpuinfo
                // map
            }
        } else { // we do subtraction according to last iteration information
            for (int t = 0; t < termIDArray.length; t++) {
                int wordID = termIDArray[t];
                topicWordCount[topic][wordID] += flag * oldWordWeight[wordID];
                sumTopicWordCount[topic] += flag * oldWordWeight[wordID];
            }
            for (int t = 0; t < termIDArray.length; t++) {
                Map<Integer, Double> gpuInfo = wordGPUInfo.get(docID).get(t);
                int wordID = termIDArray[t];
                // boolean gpuFlag = wordGPUFlag.get(docID).get(t);
                if (gpuInfo.size() > 0) {
                    for (int similarWordID : gpuInfo.keySet()) {
                        double v = gpuInfo.get(similarWordID);
                        topicWordCount[topic][similarWordID] -= v * oldWordWeight[similarWordID];
                        sumTopicWordCount[topic] -= v * oldWordWeight[similarWordID];
                    }
                }
            }
        }
    }

    public void inference() throws IOException {
//        writeDictionary();
        System.out.println("Running Gibbs sampling inference: ");

        long startTime = System.currentTimeMillis();
        System.out.println("BurnIn:");
        boolean iterFlag = false; // burnIn : iterFlag = false
        for (int iter = 1; iter <= (numIterations + numBurnIn); iter++) {
            if(iter % 50 == 0){
                if (!iterFlag) {
                    System.out.print(" " + (iter));
                } else {
                    System.out.print(" " + (iter - numBurnIn));
                }
                if (iter == numBurnIn){
                    iterFlag = true;
                    System.out.println("\nCommon iterations:");
                }
            }

            if (iter >= numBurnIn){ // Update variable when common iteration
                updateTopicProbabilityGivenWord();
                updateTopicTopWord();
                updateWordTopicSimilarity();
                if (iter % 50 == 0) {
                    updateSigma();
                    updateBDC();
                }
            }

            for (int s = 0; s < Corpus.size(); s++) {
                int[] termIDArray = Corpus.get(s);
                int preTopic = assignmentList[s];

                ratioCounter(preTopic, s, termIDArray, -1);
                for (int topic = 0; topic < numTopics; topic++) {
                    double pz = (docTopicCount[topic] + alpha) / (numDocuments - 1 + alphaSum);
                    double value = 1.0;
                    for (int t = 0; t < termIDArray.length; t++) {
                        int termID = termIDArray[t];
                        value *= (topicWordCount[topic][termID] + beta + occurenceToIndexCount.get(s).get(t) - 1)
                                / (sumTopicWordCount[topic] + betaSum + t);
                    }
                    value = pz * value;
                    multiPros[topic] = value;
                }
                int newTopic = FuncUtils.nextDiscrete(multiPros);

                // update
                assignmentList[s] = newTopic;
                ratioCounter(newTopic, s, termIDArray, +1);
            }
            if (iter >= numBurnIn && (iter % 50 == 0)) {
                updateWordWeight();
            }

        }
        expName = orgExpName;

        iterTime =System.currentTimeMillis()-startTime;
        System.out.println();
        System.out.println("Writing output from the last sample ...");

        write();
        System.out.println("Sampling completed!");

    }

    public void writeParameters() throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath + expName + ".paras"));
        writer.write("-model" + "\t" + "APUDMM");
        writer.write("\n-corpus" + "\t" + corpusPath);
        writer.write("\n-ntopics" + "\t" + numTopics);
        writer.write("\n-alpha" + "\t" + alpha);
        writer.write("\n-beta" + "\t" + beta);
        writer.write("\n-threshold" + "\t" + threshold);
        writer.write("\n-filterSize" + "\t" + filterSize);
        writer.write("\n-niters" + "\t" + numIterations);
        writer.write("\n-twords" + "\t" + top);
        writer.write("\n-name" + "\t" + expName);
        writer.write("\n-schema size" + "\t" + schemaMap.size());
        writer.write("\n-average num in schema" + '\t' + count/schemaMap.size());
        if (tAssignsFilePath.length() > 0)
            writer.write("\n-initFile" + "\t" + tAssignsFilePath);
        if (savestep > 0)
            writer.write("\n-sstep" + "\t" + savestep);

        writer.write("\n-initiation time" + "\t" + initTime);
        writer.write("\n-one iteration time" + "\t" + iterTime/(numIterations+numBurnIn));
        writer.write("\n-total time" + "\t" + (initTime+iterTime));

        writer.write("\n-word weighting : \n");
        for (int i = 0; i < vocabularySize; i++) {
            writer.write(bdcWeight[i] + " ");
        }
        System.out.println("**********************************");
        writer.close();
    }

    public void writeDictionary() throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".vocabulary"));
        for (int id = 0; id < vocabularySize; id++)
            writer.write(id2WordVocabulary.get(id) + "\n");
        writer.close();
    }

    public void writeTopTopicalWords() throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".topWords"));

        for (int tIndex = 0; tIndex < numTopics; tIndex++) {

            Map<Integer, Double> wordCount = new TreeMap<Integer, Double>();
            for (int wIndex = 0; wIndex < vocabularySize; wIndex++) {
                wordCount.put(wIndex, topicWordCount[tIndex][wIndex]);
            }
            wordCount = FuncUtils.sortByValueDescending(wordCount);

            Set<Integer> mostLikelyWords = wordCount.keySet();
            int count = 0;
            for (Integer index : mostLikelyWords) {
                if (count < top) {
                    double pro = (topicWordCount[tIndex][index] + beta)
                            / (sumTopicWordCount[tIndex] + betaSum);
                    pro = Math.round(pro * 1000000.0) / 1000000.0;
                    writer.write(id2WordVocabulary.get(index) + " ");
                    count += 1;
                }
                else {
                    writer.write("\n");
                    break;
                }
            }
        }
        writer.close();
    }

    public void writeDocTopicPros() throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".theta"));
        for (int i = 0; i < numDocuments; i++) {
            for (int tIndex = 0; tIndex < numTopics; tIndex++) {
                writer.write((pdz[i][tIndex]) + " ");
            }
            writer.write("\n");
        }
        writer.close();
    }

    public void writeTopicAssignments() throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".topicAssignments"));
        for (int dIndex = 0; dIndex < numDocuments; dIndex++) {
            int docSize = Corpus.get(dIndex).length;
            int topic = assignmentList[dIndex];
            for (int wIndex = 0; wIndex < docSize; wIndex++) {
                writer.write(topic + " ");
            }
            writer.write("\n");
        }
        writer.close();
    }

    public void writeTopicWordPros() throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(folderPath
                + expName + ".phi"));
        for (int i = 0; i < numTopics; i++) {
            for (int j = 0; j < vocabularySize; j++) {
                double pro = (topicWordCount[i][j] + beta)
                        / (sumTopicWordCount[i] + betaSum);
                writer.write(pro + " ");
            }
            writer.write("\n");
        }
        writer.close();
    }

    public void write() throws IOException {
        compute_phi();
        compute_pz();
        compute_pzd();

        writeParameters();
        writeTopTopicalWords();
        writeDocTopicPros();
//        writeTopicAssignments();
//        writeTopicWordPros();
    }

    public static void main(String[] args) throws Exception {
        String model = "KIDSDMM";

        // the parameters of model
        int topic = 50;
        double threshold = 0.7;
        double alpha = 50.0 / topic;
        double beta = 0.01;
        int lambda = 3;

        String dataName = "SearchSnippets";
        String corpusPath = "./datasets/"+ dataName + ".txt";
        String resDir = "./results/" + dataName + "/";
        String resFileName = String.format("%s_%d_%.1f_%.2f_%d_%.1f", model, topic, alpha, beta, lambda, threshold);

        KIDSDMM kidsdmm = new KIDSDMM(corpusPath, resDir, "D:\\TopicModel\\dataset\\shortProcessed\\"+dataName+"\\"+dataName+"_Word2VecSim.txt",
                threshold,20, topic, alpha, beta, lambda, 1000,
                100, 20, resFileName);
        kidsdmm.inference();
    }
}
