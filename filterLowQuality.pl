use warnings;
use strict;

if( @ARGV != 1 )
{
	die "Usage: perl filterLowQuality.pl <input> example.fastaq\n";
}
open(PTS,"$ARGV[0]") or die "Can't open file $ARGV[0]\n";
open(OUT,">highQ_$ARGV[0]");
open(LG,">>summaryStats.txt");


my($line,$line_count,@array,@sequence,$i,@quality,$num,$othername,$firstname);
$line_count=0;
while($line=<PTS>)
{
	
	$line_count ++ ;
	if($line_count%4==1)
	{
		$firstname=$line;
		@sequence=@array=@quality=();
		}
	elsif($line_count%4==3)
	{
		$othername=$line;	
		}
	elsif($line_count%4==2)
	{
		chomp($line);
		@sequence=split(//,$line);
		}	
	elsif($line_count%4==0)
	{
		chomp($line);
		@array=split(//,$line);
		for($i=0;$i<scalar(@array);$i++)
   	{
   		$quality[$i]=ord($array[$i])-33;
   	}
   	$num=scalar(@quality);
   	for($i=$num-1;$i>=0 and $quality[$i]<20;$i--)
   	{
   		$sequence[$i]='';
   		$quality[$i]='';
   		$array[$i]='';
   		}
   		if($sequence[0] eq ''){next;}
   		print OUT $firstname;	
   for($i=0;$sequence[$i];$i++)
   {
   	print OUT $sequence[$i];
   	}
   	print OUT "\n";
   	print OUT $othername;	
   for($i=0;$i<scalar(@array);$i++)	
	{
			print OUT "$array[$i]";
		}
		print OUT "\n";
	}
	}
	print LG "Raw reads	",$line_count/4,"\n";
	
close(PTS);
close(OUT);
close(LG);

