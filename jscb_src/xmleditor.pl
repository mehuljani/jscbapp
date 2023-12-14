#!/usr/bin/perl -w
#perl cgview_xml_builder.pl -sequence 2 -genes feature -title "Pseudomonas aeruginosa" -global_label T -gene_labels T -labels_to_show label -custom featureThickness=45 labelFontSize=45 -output 2_2.xml

$f=$ARGV[0];
#$f='vib.xml';
open(F,$f);
@file=<F>;
$count=0;
$f1=$ARGV[1];
#$f1='vib1.xml';
#open(F1,">$f1");
open(F1,">$f1\_2.xml");
#chmod 0777, "$f1\_1.xml";

#open(F1,">$f1\_1.xml");
foreach(@file){
	if ($_=~m/showLabel/){
		#$_=~s/showLabel="false"/label=\"GI-$count" showLabel="true"/;
		#$_=~s/showLabel="true"/label=\"GI-$count" showLabel="true"/;
		$count+=1;
	}
	
}		

foreach(@file){
	if ($_=~m/showLabel/){
		$_=~s/showLabel="false"/label=\"GI-$count" showLabel="true"/;
		
		#$_=~s/showLabel="true"/label=\"GI-$count" showLabel="true"/;
		$count-=1;
	}
	if ($_=~m/lower-center/){
		#$_=~s/lower-center/middle-center/;
	}
	if ($_=~m/plain, 20/){
		$_=~s/plain, 20/plain, 45/;
	}
	if ($_=~m/showShading="true"/){
		$_=~s/showShading="true"/showShading="false"/;
	}
	if ($_=~m/51,51,51/){
			$_=~s/51,51,51/0,0,204/;
	}
	if ($_=~m/plain, 80" text/){
			$_=~s/plain, 80" text/italics, 80" text/;
	}


	print F1 "$_";
}				

