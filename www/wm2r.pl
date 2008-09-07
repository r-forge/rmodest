#!/usb/bin/perl

$i = 0;
@model; @const; @lineA; @lineB; @a1; @b1; @c1; @s1; @a2; @b2; @c2; @s2; @msg; @value

while(<>){
	# if "Hypothesis test results: foo"
	#	if Gompertz $model[$i] = 'g' ... etc.
	# if "Lines used in this analysis: foo"
	#	$lineA[$i], $lineB[$i] = splie($_,...)
	# if "******foo"
	# 	$msg[$i] = foo;
	# if "NOTE: "
	# 	$msg[$i] .= foo;
	# if "Analysis assumes each population has unique (foo)"
	# 	if "but same (bar)"
	# 		$const[$i] = bar;
	#	else $const[$i] = 'NA';
	# if LIKELIHOOD: 
	#	$value[$i] = foo
	# if ^a: foo
	#	if ($const[$i] == 'a') $a1[$i] = $a2[$i] = foo
	#	elsif (@a1 < $i) $a1[$i] = foo
	#	else $a2[$i] == foo
	# ditto b, c, s
	# if ($model[$i] == 'g') $c1[$i] = $c2[$i] = $s1[$i] = $s2[$i] = 'NA';
	# ditto gm, l, lm
	# if "-------------------"
	# 	if (@msg < $i) $msg[$i] = 'NA';
	# 	if (the length of any array is <> $i+1, error out
	# $i++;
}

# making r vectors...
@foo = ('a','b','c');
$bar = join(",",@foo);
print("foo<-c($bar)\n");