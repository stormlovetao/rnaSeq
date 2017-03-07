
#!/bin/bash

for x in /PHShome/tw786/neurogen/Tao/kraken_standard_db/library/Viruses/*; do
	count=`ls -l $x/*.fna 2> /dev/null | wc -l`
	if [[ $count == 0 ]]; then
		echo $x
		cat $x/*.ffn >> ~/virus.fa
	else
		cat $x/*.fna >> ~/virus.fa
	fi
done