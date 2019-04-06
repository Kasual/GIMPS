<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<!-- $Id$ -->
<html>

<head>
  <title>Glucas -Yet Another FFT- spanish</title>
  <meta name="GENERATOR" content="Quanta Plus">
  <meta name="AUTHOR" content="Guillermo Ballester Valor">
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-15">
  <meta name="KEYWORDS" content="FFT, convoluci�n, Transformada R�pida de Fourier, N�meros de Mersenne,
   test de Lucas-Lehmer">
  <meta name="DESCRIPTION" content="Librer�a FFT para hacer convoluciones tan r�pido como sea posible. Glucas es 
  una implementaci�n del test de Lucas-Lehmer para determinar la primalidad de los n�meros de Mersenne. Utiliza
  la librer�a YEAFFT">
</head>
<body lang="es-ES" text="#000000" link="#0000ef" vlink="#51188e"
 bgcolor="#ccffff">
<center>
<img src="http://sourceforge.net/sflogo.php?group_id=24518&amp;type=5"
 name="Imagen1" alt="SourceForgeLogo" align="bottom" width="210"
 height="62" border="0">

<a href="http://petition.eurolinux.org"><img src="../patent_banner.png"></a>
</center>
 
 
<h1><b><i><font size="7"><font face="Arial, Helvetica"><font
 color="#993300">G</font></font></font></i>LUCAS - Yet Another FFT</b></h1>
<h3><a href="./index.html">Introducci�n en Ingl�s</a>  </h3>

<h3><a href="#Introducci�n">Introducci�n</a> </h3>� � � � � � � � � � � � � �
 
<h3><a href="http://sourceforge.net/projects/glucas">Estatus del proyecto (sourceforge)</a> </h3>
 
<h3><a
 href="http://sourceforge.net/project/showfiles.php?group_id=24518">Bajar (vieja version  2.9.0)</a> </h3>
 
<h3><a
 href="../pub/glucas/snapshots">�ltimas versiones tar.gz (version 2.9.2)</a> </h3>

<h3><a href="http://www.oxixares.com/viewcvs">Navegar por el repositorio</a></h3>

<h4><a href="http://glucas.sourceforge.net/oldbin.html">Binarios a�n m�s antiguos</a></h4>
 
<h3><a href="./doc/index.html">Documentaci�n en l�nea (Ingl�s)</a></h3>
 

<p><a name="Introducci�n"></a>I<b>NTRODUCCI�N</b>  </p>
 
<p><font face="Arial, Helvetica"><b>Glucas </b>es un programa de c�digo abierto
y gratuito para realizar el test primalidad de los n�meros de Mersenne (n�meros
de la forma 2^n - 1). Sus resultados nos dicen si un n�mero de Mersenne es
un n�mero primo. Puede leer m�s acerca de los n�meros de Mersenne <a
 href="http://www.utm.edu/research/primes/mersenne.shtml">aqu�.</a> La forma
especial de los n�meros de Mersenne ofrece tres grandes ventajas:</font>
 </p>
 
<p><font face="Arial, Helvetica">1) Existe un test suficientemente r�pido
(test de Lucas-Lehmer). El resultado de dicho test es positivo si y solo
si el n�mero de Mersenne es primo.</font>  </p>
 
<p><font face="Arial, Helvetica">2) La aritm�tica modular que se precisa
en los c�lculos puede realizarse de una forma sencilla sin tener que recurrir
a costosas divisiones.</font>  </p>
 
<p><font face="Arial, Helvetica">3) El algoritmo utilizado para multiplicar
los n�meros tan grandes que necesitamos (algunos millones de bits) utilizan
una forma especial de Transformada Rapida de Fourier, <a
 href="http://perfsci.com/">DWT</a> . La velocidad que se alcanza con esta
Transformada es el doble de la forma general de multiplicar grandes enteros
a traves de Transformadas R�pidas de Fourier y el teorema de convoluci�n.</font>
 </p>
 
<p><font face="Arial, Helvetica">Estas ventajas hacen que la l<a
 href="http://www.utm.edu/research/primes/largest.html">ista de records de
grandes numeros primos</a> est� liderada por los n�meros de Mersenne. Usted
puede ser el descubridor de un nuevo n�mero primo de Mersenne record, pero
no se haga ilusiones r�pidamente, es bastante dif�cil. Hay un proyecto en
Internet dedicado a la b�squeda de esos n�meros tan escasos. Este proyecto
es <a href="http://mersenne.org/prime.htm">GIMPS, the Great Internet Mersenne
Prime Search</a></font> . <font face="Arial, Helvetica">GIMPS ya ha descubierto
los cuatro n�meros primos m�s grandes de la historia y es un proyecto pionero
en el c�lculo distribuido por Internet. El n�cleo del proyecto es G.Woltman, 
quien ha escrito los programas clientes para la plataforma X86. Sus programas
son desconcertantemente r�pidos, y pienso que es la transformada r�pida de
Fourier m�s r�pida que se ha escrito para dicha plataforma. Esta escrita
cuidadosamente optimizada en lenguaje ensamblador.</font>  </p>
 
<p><font face="Arial, Helvetica">Glucas es, en cambio, un programa escrito
en C dise�ado para realizar el test de Lucas-Lehmer. Como yo todav�a trabajo
con antiguos pentiums en casa y en el trabajo, tambi�n he incluido algunas
macros en ensamblador para esta plataforma utilizando el excelente compilador
GNU/gcc. Se ha alcanzado un 80% de la velocidad de la �ltima version publicada
de G.Woltman (prime95 v.20). Para otras plataformas solamante hay dos clientes
suficientemente r�pidos, Glucas y Mlucas. <a
 href="http://novarese.net/prime/mlucas.html">Mlucas</a> es un excelente
programa escrito en fortran-90 escrito por Ernst W. Mayer (una adaptaci�n
de Mlucas se utiliza en el test specf2000). La velocidad de Mlucas es similar
a Glucas. Mlucas necesita un compilador f-90 y Glucas un compilador C. Para
algunas plataformas como MAC, no hay buenos compiladores fortran-90 disponibles
y aqu� radica la ventaja de Glucas.</font>  </p>
 
<p><font face="Arial, Helvetica">Puede leer m�s acerca de Glucas. Tiene a
su disposici�n el manual completo en l�nea <a
 href="http://glucas.sourceforge.net/glucas/index.html">aqu�</a>. Tambi�n
puede obtener el <a
 href="http://prdownloads.sourceforge.net/glucas/Glucas-2.9.0.docs.tar.gz?download">paquete
de documentaci�n</a> que incluye el manual en formatos HTML, DVI y PDF.</font></p>
 
<p><font face="Arial, Helvetica"><b>YEAFFT</b> es la librer�a utilizada por
Glucas para las convoluciones (grandes multiplicaciones) que se necesitan.
Esta escrita por m�, Guillermo Ballester Valor. Est� basado <a
 href="http://www.perfsci.com/free/techpapers/confgt.pdf">este art�culo</a>
. Es una librer�a escrita en C con algunas macros definidas en ensamblador
. Se ha escrito bajo licencia GPL.</font> </p>
 
<p><br>
<br>
<h3>ACTUALIZACI�N: Debido a algunas dificultades de acceso en Sourceforge, estoy moviendo parcialmente 
el proyecto a mi propio servidor. La versi�n actual en desarrollo es la 2.9.2. Puede bajarse
lo �ltimo <a href="../pub/glucas/snapshots"> aqu� </a>
<br>
</p>
<p>Ademas del repositorio CVS en Sourceforge, puede acceder de forma an�nima al 
repositorio tipo <a href="http://subversion.tigris.org">Subversion</a>. Para bajarse el repositorio:</p>
<CODE>svn co https://svn.oxixares.com/repos/glucas</CODE>
<p>Para actualizarse con las �ltimas modificaciones:</p>
<CODE>svn update https://svn.oxixares.com/repos/glucas</CODE>
<p>
El drectorio actualmente activo en el repositorio es 'glucas/trunk/glucas'
</p>
<p>Tambi�n puede <a href="https://svn.oxixares.com/viewcvs">ver el c�digo fuente con su navegador</a></p>
<code>https://svn.oxixares.com/viewcvs</code>


 <br>
<h5><b>ACTUALIZADO: 06-Mar-2005</b></h5>

</body>
</html>
