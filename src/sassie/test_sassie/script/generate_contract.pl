#!/usr/bin/perl
sub ltrim($)
{
   my $string = shift;
   $string =~ s/^\s+//;
   return $string;
}

$file = $ARGV[$0];
open(FILE,$file);
@lines=<FILE>;
$t=0;
foreach $line (@lines)
{
   if ($line =~ "'''")
   {
      if ($t==0) {$t=1;}
      else
      {
         if($t==1) {$t=0;}
      }
   }
   else
   {
      if ($t==1) {print ltrim($line);}
   }
}
