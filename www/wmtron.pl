use File::Copy;
$mini = 0;
$maxi = 100;
$startj = 1;
$maxj = $maxi-1;
$precision = "10";
$checkboxht = "455";
$wmpath = "E:\\WinModest\\Exampl~1\\";
$wmcmd = "E:\\WinModest\\WinModest.exe";
$wmdat = "sim100.dat";
$wmout = "sim100.out";
$begintime = time;
$numrecords = 0;

@down1 = (0,1,2,3);
@tab1n2 = (2,3,3,4);
@tab3 = (9,11,11,13);
@tab4 = (6,7,7,8);

#@badlines = (57, 58, 65, 80, 83, 93, 97, 98);

@starts = (["1e-6","1e-6","1e-6","1e-6"],
       ["1e-6","1e-6","1e-6","1e-6","1e-6","1e-6"],
       ["1e-6","1e-6","1e-6","1e-6","1e-6","1e-6"],
       ["1e-6","1e-6","1e-6","1e-6","1e-6","1e-6","1e-6","1e-6"]);

@models = ([455,495],[455,495,565],[455,495,535],[455,495,565,535]);

if (open(INPUT,"e:\\winmodest\\wmcrash")){
    @crashrec = <INPUT>;
    $mini = $crashrec[0]; chomp($starti);
    $startj = $crashrec[1]; chomp($startj);
    print "File contents: $starti , $startj\n";
}
close(INPUT);

#print("starting\n");
#exec($wmcmd); print("started\n");
#system("tron type {wait}{alt}foa$wmdat{enter}{wait}");
#system("tron type {alt}{tab}");

for $i ($mini .. $maxi){
    print "comparing $i to: ";
    if($mini>0 & $i == $mini){$z = $startj;}else{$z = 1;}
    for $j ($i+$z .. $maxj){
        #($imatch) = grep /$i/ ,@badlines;
        #($jmatch) = grep /$j/ ,@badlines;
        #if($imatch|$jmatch){print("Bad population/s ($i vs $j).\n");next;}else{
        print "$j ... ";
        for $k (0 .. $#models){
            $model = $models[$k];
            $start = $starts[$k];
            $par2start = $tab1n2[$k];
            $par2end = 2*$par2start-1;
            $par1end = $par2start -1;
            $n = @$model -1;
            for $l (0 .. $n){
                $win = `tron GetWindow`; chomp($win);
                if( $win !~ /WinModest/ ){ 
                    open $crash, '>', "e:\\winmodest\\wmcrash" or die "Can't write file";
                    print $crash "$i\n";
                    print $crash "$j\n";
                    close $crash;
                    copy($wmpath.$wmout, $wmpath.$wmout.$i.vs.$j);
                    die "\nWinModest crashed comparing $i to $j. Restarting.\n";
                }
                $xcoord = $$model[$l];
                $cmd = "tron type {alt}mmh{tab}{tab}";
                $cmd .= ("{down}") x $down1[$k];
                for $m (0 .. $par1end){
                    $cmd .= "{tab}$$start[$m]";
                }
                $cmd .= "{tab}{end}{tab}".$precision."{tab}{tab}";
                for $m ($par2start .. $par2end){
                    $cmd .= "{tab}$$start[$m]";
                }
                $cmd .= "{tab}{end}{enter}";
                $cmd .= ("{tab}") x ($tab3[$k]);
                $cmd .= "{home}";
                $cmd .= ("{down}") x ($i);
                $cmd .= ("{tab}") x ($tab4[$k]);
                $cmd .= "{home}";
                $cmd .= ("{down}") x ($i);
                $cmd .= "{enter}";
                $cmd .= ("{down}") x ($j-$i-1);
                system($cmd);
                #if ($xcoord > 0) { system("tron leftclick $xcoord,$checkboxht"); }
                system("tron leftclick $xcoord,$checkboxht");
                system("tron type {enter}{wait}{wait}");
            }
        }
        system("tron type {alt}fw{enter}{wait}y");
        $numrecords++; $runtime = time - $begintime;
        print "\n $numrecords total comparisons made in $runtime\n";
        print "Average speed: ".$runtime/$numrecords." seconds per record.\n";
    }
    copy($wmpath.$wmout, $wmpath.$wmout.$i.vs.$j) or print "Copy failed: $!";
}

$numrecords++; $runtime = time - $begintime;
open $done, '>', "e:\\winmodest\\wmdone" or die "Can't write file wmdone";
print $done, "WinModest run successfully completed.\n";
print $done, "$numrecords total comparisons made in $runtime\n";
print $done, "Average speed: ".$runtime/$numrecords." seconds per record.\n";
close $done;
