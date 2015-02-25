#Set CLASSPATH to contemplate.jar and dependencies
export CLASSPATH=../lib/contemplate.jar:../lib/nmica.jar:../lib/biojava.jar:../lib/bytecode.jar

#Make submodel XML files
create:
	java MakeBlock 2 .5 .5 .1 > block.xml
	java MakeType1spacer 1 1 .5 .5 .1 > type1spacer.xml
	java MakeType2spacer 1 1 .5 .5 .1 > type2spacer.xml

#Train submodels using Viterbi algorithm
train:
	java TrainBlock block.xml blocks.fa 2 > block_train.xml
	java TrainType1spacer type1spacer.xml type1spacers.fa 2 > type1spacer_train.xml
	java TrainType2spacer type2spacer.xml type2spacers.fa 2 > type2spacer_train.xml
	##Construct full Contemplate enhancer search model
	java MakeEnhSearch block_train.xml type1spacer_train.xml type2spacer_train.xml > enhSearch.xml

#Run enhancer search model to detect CRMs
test:
	java -server -Xmx1300M ModelRunner enhSearch.xml tll_mel.fa tll_vir.fa -fastaOut tll.mfa -blockListOut tll.blk -laganOut tll.lgn
