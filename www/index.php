
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<! --- R-Forge Logo --- >
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p>Breaking news. Survomatic is now available as a standard package. To install, type: install.packages("survomatic", repos = "http://r-forge.r-project.org")<br>
But be warned, this is my first package ever, and before I packaged it, I was distributing Survomatic using a system of my devising, and there may be some legacy assumptions the scripts make that are no longer true. These will be fixed quickly, because they should for the most part be glaringly obvious.</p>

<ul>
<li><a href='rmodest3.r'>rmodest3.r</a>, the core mortality estimation module.</li>
<li><a href='survomatic.r'>survomatic.r</a>, a user friendly wrapper for rmodest3 and a whole raft of other survival functions. Its default behavior at the moment is to prompt you for everything and then output a self-executing .rdata file that you can run again later. Produces Kaplan-Meier survival plots, quantile regression plots, and mortality hazard plots.</li>
<li><a href='simsurv.r'>simsurv.r</a>, a script for creating arbitrary numbers of datasets with known Gompertz/Logistic parameters. Useful for testing survival analysis software.</li>
</ul>

<p>More to come, stay tuned.</p>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
