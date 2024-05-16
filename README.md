# KIDS-DMM (Keyword Informative Discriminative Scheme Dirichlet Multinomial Mixture)

KIDS-DMM uses an [open-source Java package](https://github.com/qiang2100/STTM) to implement the algorithm.

## 1. Requirements

- Java （Version=1.8）

## 2. Datasets
We procided the following three short text datasets for evaluation, SearchSnippets, GoogleNews, and Biomedical. All of corpus files and the corresponding label files have been prepared in the path **./datasets** according to the survey, [Short Text Topic Modeling Techniques, Applications, and Performance: A Survey](https://ieeexplore.ieee.org/document/9086136). 

Taking SearchSnippets as an example, the dataset file path is as follows.

>datasets
>> SearchSnippets
>>> word_wiki
>>> 
>>> SearchSnippets.txt
>>> 
>>> SearchSnippets_label.txt
>>> 
>>> SearchSnippets_vocab.txt
>>>
>>> SearchSnippets_Word2VecSim.txt

For the corresponding word_wiki and the word2VecSim, you can download from [this](https://drive.google.com/drive/folders/1Vw1OOedBjP5c93HJIZZGbCy100q6kOgn) following this [paper](https://link.springer.com/chapter/10.1007/978-981-99-8181-6_28).

## 3. Run and Evaluate KIDS-DMM
    bash run.sh