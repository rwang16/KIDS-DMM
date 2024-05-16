current_dir=$(pwd)

project_jar=${current_dir}/out/artifacts/KIDS-DMM.jar

algorithm=KIDSDMM
K=50
alpha=1.0
beta=0.01
lambda=3

# corpus_name:{SearchSnippets, GoogleNews, StackOverflow}
# the corresponding threshold :{0.7, 0.8, 0.8}
corpus_name=SearchSnippets
threshold=0.7

corpus_path=${current_dir}/datasets/${corpus_name}/${corpus_name}.txt
word_sim_path=${current_dir}/datasets/${corpus_name}/${corpus_name}_Word2VecSim.txt
res_dir=${current_dir}/results/${corpus_name}/

res_file_name=${algorithm}_${K}_${alpha}_${beta}_${threshold}_${lambda}

java -jar ${project_jar}\
     -model ${algorithm}\
     -dataname ${corpus_name}\
     -alpha ${alpha}\
     -beta ${beta}\
     -lambda ${lambda}\
     -ntopics ${K}\
     -schema ${word_sim_path}\
     -threshold ${threshold}\
     -corpus ${corpus_path}\
     -output ${res_dir}\
     -name ${res_file_name}\
     -nBurnIn 100\
     -niters 1000

# Topic quality evaluation
vocab_path=${current_dir}/datasets/${corpus_name}/${corpus_name}_vocab.txt
wiki_dir=${current_dir}/datasets/${corpus_name}/word_wiki
tc_evaluation=pmi

java -jar ${project_jar}\
      -model "TopicQualityEval"\
      -ntopics ${K}\
      -topWordsDir ${res_dir}\
      -vocabFile ${vocab_path}\
      -wikiDir ${wiki_dir}\
      -topicCoherEval ${tc_evaluation}\
      -topTC 15

# Classification evaluation
label_path=${current_dir}/datasets/${corpus_name}/${corpus_name}_label.txt
java -jar ${project_jar}\
      -model ClassificationEval\
      -label ${label_path}\
      -dir ${res_dir}\
      -prob theta