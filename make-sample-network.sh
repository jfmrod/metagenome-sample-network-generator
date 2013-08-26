#!/bin/bash

OTHERSAMPLES=""
SAMPLESATTR=""

while getopts "o:a:h" OPT; do
  case $OPT in
    o) OTHERSAMPLES=$OPTARG ;;
    a) SAMPLESATTR=$OPTARG ;;
    h) cat <<EOF
make-sample-network.sh v1.0

Script to generate a network of metagenomic samples
Network is produced in dot file format to be visualized with neato (graphviz)

EOF
exit -1
;;
  esac
done

shift $(( OPTIND-1 ))

if [ -z "$1" ]; then
  cat <<EOF
make-sample-network.sh v1.0

Script to generate a network of metagenomic samples
Network is produced in dot file format to be visualized with neato (graphviz)

syntax:
  make-sample-network.sh [-a samples.attrs] [-o othersamples.otu] <sample1.otu> [sample2.otu] [...]

help:
  make-sample-network.sh -h
EOF
  exit -1
fi

SAMPLES=$*


awk '
BEGIN{
  othercutoff=0.20;
  othercutoff2=0.10;
  samplecutoff=0.25;

  FS="\t";
  f="'$OTHERSAMPLES'";
  if (length(f)>0){
    while (getline<f > 0){
      if (/^#/) continue;
      if (/^>/){ otu=substr($1,2); otus[otu]=1; continue; }
      accotu[$1]=otu;
      if (!(otu SUBSEP $2 in otusample)){
        sampleotus[$2]=sampleotus[$2] SUBSEP otu;
        ++sampleotucount[$2];
      }
      ++otusample[otu SUBSEP $2];
      othersampletype[$2]=1;
    }
  }
  f="'$SAMPLESATTR'";
  if (length(f)>0){
    while (getline<f > 0){
      if (/^#/) continue;
      if (/^>/) { sid=substr($1,2); continue; }
      split($0,a,"=");
      samplesinfo[sid SUBSEP a[1]]=a[2];
      samplesattr[a[1]]=1;
    }
  }
  FS=" ";
  print "graph G {";
  print "size=\"6,3!\";";
  print "node [fontname=\"Arial\",fontsize=\"32.0\",width=0.5,penwidth=0];";
  print "outputorder=edgesfirst;";
}

function min(a,b){return(a<b?a:b);}

function getBasename(f,      a)
{
  split(f,a,"/");
  return(a[length(a)]);
}

/^#/ { next; }

{
  type=getBasename(FILENAME);
  sampletype[type]=1;
  otu=$2;
  otus[otu]=1;

  if (!(otu SUBSEP type in otusample)){
    sampleotus[type]=sampleotus[type] SUBSEP otu;
    ++sampleotucount[type];
  }
  ++otusample[otu SUBSEP type];
}

END{
  print "# other samples: " length(othersampletype);
  for (i in sampletype)
    print "# sample: " i " otus: " sampleotucount[i];
  for (i in othersampletype){
#    if (sampleotucount[i]<20) continue;
    ++totalotherstudies;
    print "# other sample: " i " otus: " sampleotucount[i];
  }
  print "# total other samples (>=20 OTUs): " totalotherstudies;

  for (i in sampletype){
    for (j in othersampletype){
#      if (sampleotucount[j]<20) continue;
      split(sampleotus[j],a,SUBSEP);
      for (k=2; k<=length(a); ++k){
        if ((a[k] SUBSEP i) in otusample) ++sharedotu[i SUBSEP j];
      }
      f=sharedotu[i SUBSEP j]/min(sampleotucount[i],sampleotucount[j]);
      if (f>=othercutoff2)
        otherotu[j]=1;
    }
  }

  samplecount=0;
  for (i in sampletype){
    delete attr;
    attr["shape"]="box";
    attr["color"]="\"#FF3333\"";
    attr["fillcolor"]="\"#FF3333\"";
    attr["label"]="\"" i "\"";
    attr["width"]=0.6;
    attr["height"]=attr["width"];
    attr["fixedsize"]="true";
    attr["style"]="filled";
    for (j in samplesattr){
      if (i SUBSEP j in samplesinfo)
        attr[j]=samplesinfo[i SUBSEP j];
    }
    attrstr="";
    for (j in attr)
      attrstr=attrstr "," j "=" attr[j];

    print "n" samplecount " [" substr(attrstr,2) "];";
    sampleid[i]=samplecount;
    ++samplecount;
  }
  for (i in otherotu){
    delete attr;
    attr["shape"]="circle";
    attr["color"]="\"#000000\"";
    attr["fillcolor"]="\"#FF9999\"";
    attr["label"]="\"\"";
    attr["width"]=0.4;
    attr["height"]=attr["width"];
    attr["penwidth"]=0;
    attr["fixedsize"]="true";
    attr["style"]="filled";
    for (j in samplesattr){
      if (i SUBSEP j in samplesinfo)
        attr[j]=samplesinfo[i SUBSEP j];
    }
    attrstr="";
    for (j in attr)
      attrstr=attrstr "," j "=" attr[j];

    print "n" samplecount " [" substr(attrstr,2) "];";
    sampleid[i]=samplecount;
    ++samplecount;
  }
  for (i in sampletype){
    for (j in otherotu){
      f=sharedotu[i SUBSEP j]/min(sampleotucount[i],sampleotucount[j]);
      if (f>othercutoff)
        printf "n" sampleid[i] " -- n" sampleid[j] " [len=%f,weight=%f,penwidth=%f,color=\"#000000%2x\"];\n",0.3*(1.0-f)+4.0,10.0*(1.0-f),(f-othercutoff)*3.0/(1.0-othercutoff)+1.5,int(((f-othercutoff)*0.30/(1.0-othercutoff)+0.1)*255);
    }
  }
  for (i in sampletype){
    for (j in sampletype){
      if (i==j) break;
      for (k in otus){
        if ((k SUBSEP i) in otusample && (k SUBSEP j) in otusample) ++sharedotu[i SUBSEP j];
        if ((k SUBSEP i) in otusample || (k SUBSEP j) in otusample) ++totalotu[i SUBSEP j];
      }
      f=sharedotu[i SUBSEP j]/min(sampleotucount[i],sampleotucount[j]);
      if (f>samplecutoff){
        tmp=int((0.7+(f-samplecutoff)*0.3/(1.0-samplecutoff))*255);
        printf "n" sampleid[i] " -- n" sampleid[j] " [len=%f,weight=%f,penwidth=%f,color=\"#00FF00%2x\"];\n",1.0*(1.0-f)+3.0,20.0*(1.0-f),(f-samplecutoff)*3.0/(1.0-samplecutoff)+1.5,tmp;
      }
    }
  }
  print "}";
}
' $SAMPLES

