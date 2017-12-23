#!/usr/bin/perl

# #############################################################################
# USAGE:
# perl Helitron-groups_difinition_V6.pl blast_out clases.tsv blastOUT.tsv
# #############################################################################

####################################################################################################
#*****************************
#                            	Subroutines
#                                            *******************************************************
####################################################################################################

sub search_elem{
    $Flag = 0;
    local($Start1, $End1, $Start2, $End2);
    $Start1 = $_[0];  # inicio del elemento que busca compañia
    $End1   = $_[1];  # final del elemento que busca compañia
    $Start2 = $_[2];
    $End2   = $_[3];

    $Length1 = $End1 - $Start1 + 1;
    $Length2 = $End2 - $Start2 + 1;

    if($Length1 > $Length2){
	$value      = $Length1 / $Length2;
	$length80   = $Length1 * $overlapVal;
	$lengthAway = $Length1 * $outDistance;
    } else{
	$value      = $Length2 / $Length1;
	$length80   = $Length2 * $overlapVal;
	$lengthAway = $Length2 * $outDistance;
    }

    # val <= covVal... analiza 250/100 y 250/220, mientras más pequeño sea el framento chico..
    # más grande tenderá a ser val... lo ideal es tener valores de val que oscilen en 1.. =)

    if($Start2 <= $Start1 && $End2 >= $End1){
	if(($End1 - $Start1) >= $length80 && $value <= $covValue){
	    $Flag=1;
	    # nos dice que el segmento está contenido dentro de otro segmento...
	}
    } elsif($Start1 <= $Start2 && $End1 >= $End2){
	if(($End2-$Start2) >= $length80 && $value <= $covValue){
	    $Flag=2;
	    # nos dice que el segmento está contenido dentro de otro segmento...
	}
    } elsif($Start1 <= $Start2 && $End1 > $Start2 && $End1 <= $End2){
	if(($End1-$Start2) >= $length80 && $value <= $covValue && ($Start2-$Start1) <= $lengthAway
	   && ($End2-$End1) <= $lengthAway){
	    $Flag=3;
	    # Si el sector va acá.. debemos de ser capaces de mover los limites del nuevo grupo..
	}
    } elsif($Start2 <= $Start1 && $End2 > $Start1 && $End2 <= $End1){
	if(($End2-$Start1) >= $length80 && $value <= $covValue && ($Start1-$Start2) <= $lengthAway
	   && ($End1-$End2) <= $lengthAway){
	    $Flag=4;
	    # Si el sector va acá.. debemos de ser capaces de mover los límites del nuevo grupo..
	}
    }

    # regresamos la misma variable. 1 indica que el segmento cumple con las reglas establecidas
    # 0 indica que el segmento no cumple con las reglas establecidas
    return($Flag);
}


# subrutina para la busqueda de que todos los helitrones de un grupo cumplan con el mismo requisito

sub control_param{
    local($valS, $valE, $Ky1, $Ky2);
    $valS = $_[0];
    $valE = $_[1];
    $Ky1  = $_[2];
    $Ky2  = $_[3];

    # por cada TE que conforma la lista del segmento, debemos de preguntar si se respentan las
    # mismas reglas de "conservación" entre los helitrones y el nuevo elemento identificado

    @data=split(/-/, $listElement{$Ky1}{$Ky2});  # cada uno de los helitrones que conforman el grupo
    $Flag2 = 3;                                  # bandera que nos dirá si no se cumplen los
    # requisitos de filtro con algún elemento..
    $j = 0;
    # recorremos cada uno de los helitrones.. y tomamos sus coordenadas para comparar
    {
	do{
	    $data[$j]=~ /(AT\dTE\d+)\.(\d)/;
	    $indeX = $2;
	    $localTE = $1;
	    $segment{$Ky1}{$localTE}{$indeX}=~ /(\d+)\-(\d+)/;
	    $localS = $1;
	    $localE = $2;
	    $j++;
	    $Flag2 = &search_elem($valS,$valE,$localS,$localE);   # podemos mandar la información de
	                                                           # de las coordenadas de inicio-fin
	    last if(!$Flag2);
	    # resultados posibles:
	    # 0 = no cumple, 1 = dentro, 2 = dentro, 3 = por un lado, 4 = por otro lado
	} while($j < scalar(@data) || !$Flag2);
    }  # agregamos este corchete por fuera.. para tener el loop y poder romperlo con un "last"

    return($Flag2);
}

####################################################################################################
#***************************************
#                                            Main program
#                                                                ***********************************

#---------------------------------------------------------------------------------------------------
# El propósito del script es identificar las capturas génicas en el conjunto de helitrones que
# poseen las características estructurales.
# Toma como archivo de entrada el resultado del BLAST. Analiza los alineamientos en búsqueda de
# aquellos que posean los valores y criterios definidos.
# Archivo de salida es una tabla con la información de las agrupaciones encontradas.

open(BLAST, "$ARGV[0]") or die "error i could not load the first file\n"; # el resultado del blast

$overlapVal  = 0.8;   # grado de sobrelape aceptado para considerar que dos regiones son la misma
$outDistance = 0.2;   # distancia aceptada de no sobrelape, para seguir considerando que dos
                      # regiones son las mismas, del elemento más largo??
$covValue = 1.2;      # valor de cobertura. de un segmento pequeño dentro de un segmento más grande

####################################################################################################
# OPENING BLAST OUTPUT

while(<BLAST>){
    chomp $_;
    @data = split(/\t/,$_);
    # solo aceptamos regiones de match con una longitud mayor o igual a 50bp
    next if($data[3] <= 50);
    $data[1] =~ /(\w+)\|(\+|\-)/;
    $tsense  = $2;
    $data[1] =$1;

    $data[0] =~ /(\w+)\|(\+|\-)/;
    $data[0] = $1;
    $qsense  = $2;

    # Establecemos el orden 5'-3', la coordenada de inicio es más pequeña que la del fin.
    # Esto para el target (captura génica)
    if($data[8] < $data[9]){
	$startt  = $data[8];
	$endt    = $data[9];
	$ttsense = "+";
    } else{
	$startt  = $data[9];
	$endt    = $data[8];
	$ttsense = "-";
    }
    # Establecemos el orden 5'-3' para el query (transposon)
    if($data[6] < $data[7]){
	$startq  = $data[6];
	$endq    = $data[7];
	$qtsense = "+";
    } else{
	$startq  = $data[7];
	$endq    = $data[6];
	$qtsense = "-";
    }

    ################################################################################################
    # Tenemos que soluccionar el mejor match de captura...
    # generamos arreglos donde la primer llave es el identificador del helitron y la segunda es un
    # indice .. que simplemente cuenta cuantas capturas va teniendo cada helitron...

    if(!exists($helitronIndex{$data[0]})){
	$helitronIndex{$data[0]} = 1;
	$index = $helitronIndex{$data[0]};
	# creamos HelitronEval.. para la selección de la mejor región local
	$helitronEval{$data[0]}{$index} = $data[10];

	# creamos un arreglo 'fantasma', con la información "original" y natural del blast output.
	$helitronMapGhost{$data[0]}{$index} = $data[1]."|".$data[3]."|".$data[6]."-".$data[7]."|".$data[8]."-".$data[9]."|".$qsense."|".$tsense."|".$qtsense."|".$ttsense;

	$helitronMap{$data[0]}{$index} = $data[1]."|".$data[3]."|".$startq."-".$endq."|".$startt."-".$endt;
	# hash{helitron}{indice} = genmapeo.|.longitud.|.inicio-fin helitron.|.inicio-fin genmapeo

    } else {
	# No es un helitron nuevo el que estamos observando..
	# Preguntamos si la región que estamos analizando sobrelapa con alguna región que habíamos
	# observado previamente
	# Consideramos el eval del alineamiento como factor de discriminación...

	$flag = 0;
	for $key(sort {$a <=> $b} keys %{$helitronMap{$data[0]}}){
	    @seq      = split(/\|/,$helitronMap{$data[0]}{$key});
	    $seq[2]   =~ /(\d+)\-(\d+)/;
	    $startseq = $1;
	    $endseq   = $2;

	    # el segmento que ya teníamos es mas largo. ahora preguntamos por el eval, si tiene un
	    # mejor eval cambiamos y eliminamos al anterior
	    if($startseq <= $startq && $endseq >= $endq){
		$flag = 1;
		if($data[10] < $helitronEval{$data[0]}{$key}){
		    $helitronMap{$data[0]}{$key} = $data[1]."|".$data[3]."|".$startq."-".$endq."|".$startt."-".$endt;
		    $helitronMapGhost{$data[0]}{$key} = $data[1]."|".$data[3]."|".$data[6]."-".$data[7]."|".$data[8]."-".$data[9]."|".$qsense."|".$tsense;
		    $helitronEval{$data[0]}{$key} = $data[10];
		}
	    } elsif($startq < $startseq && $endq > $endseq){
		# el segmento nuevo es más largo. consideramos el eval, cambiamos el flag quiere
		# decir que evaluamos pero que no fue funcional
		$flag=1;
		if($data[10] < $helitronEval{$data[0]}{$key}){
		    # el eval es menor.. cambiamos de sector
		    $helitronMap{$data[0]}{$key} = $data[1]."|".$data[3]."|".$startq."-".$endq."|".$startt."-".$endt;
		    $helitronMapGhost{$data[0]}{$key} = $data[1]."|".$data[3]."|".$data[6]."-".$data[7]."|".$data[8]."-".$data[9]."|".$qsense."|".$tsense;
		    $helitronEval{$data[0]}{$key} = $data[10];
		}
	    } elsif(($startseq <= $startq && $startq < $endseq && $endseq < $endq) || ($startq < $startseq && $endq > $startseq && $endq <= $endseq)){
		# ahora preguntamos si existe sobrelape...
		$flag = 1;
		if($data[10] < $helitronEval{$data[0]}{$key}){
		    # el eval es menor.. cambiamos de región por esta.. =)
		    $helitronMap{$data[0]}{$key} = $data[1]."|".$data[3]."|".$startq."-".$endq."|".$startt."-".$endt;
		    $helitronMapGhost{$data[0]}{$key} = $data[1]."|".$data[3]."|".$data[6]."-".$data[7]."|".$data[8]."-".$data[9]."|".$qsense."|".$tsense;
		    $helitronEval{$data[0]}{$key} = $data[10];
		}
	    }
	}
	if(!$flag){
	    # el elemento no muestra sobrelape.. consideramos al segmento como una nueva captura...
	    # definimos una nueva región..
	    $helitronIndex{$data[0]}++;
	    $index = $helitronIndex{$data[0]};
	    $helitronEval{$data[0]}{$index} = $data[10];
	    $helitronMapGhost{$data[0]}{$index} = $data[1]."|".$data[3]."|".$data[6]."-".$data[7]."|".$data[8]."-".$data[9]."|".$qsense."|".$tsense."|".$qtsense."|".$ttsense;
	    $helitronMap{$data[0]}{$index} = $data[1]."|".$data[3]."|".$startq."-".$endq."|".$startt."-".$endt;
	    # hash{helitron}{indice} = genmapeo.|.longitud.|.inicio-fin helitron.|.inicio-fin genmapeo
	}
    }
}

close(BLAST);

####################################################################################################
# BUSQUEDA DEL MEJOR ALINEAMIENTO HELITRON-CAPTURA GÉNICA
# creemos que en ocasiones existen ciertas regiones de homología entre genes que podrían dar
# multiples alineamientos posibles pero solo una es la mejor.
# En base de buscar al alineamiento más largo/mejor, evitamos esta redundancia y ruido a la hora de
# formar los grupos de ancestría.

# el mejor alineamiento es el mas largo.
# solo nos enfocamos dentro de las coordenadas de la captura génica dentro del helitron

undef %ehelitronIndex;
for $key(sort keys %helitronMap){

  REGRESO:
    for $key2(sort keys %{$helitronMap{$key}}){  # comenzamos con el hit del gen 1
	$helitronMap{$key}{$key2} =~ /(\w+)\|(\d+)\|(\d+\-\d+)\|(\d+\-\d+)/;
	$query  = $3;
	$query  =~ /(\d+)\-(\d+)/;
	$start1 = $1;
	$end1   = $2;

	for $key3(sort keys %{$helitronMap{$key}}){
	    if($key3 > $key2){
		$helitronMap{$key}{$key3}=~ /(\w+)\|(\d+)\|(\d+\-\d+)\|(\d+\-\d+)/;
		$query = $3;
		$query =~ /(\d+)\-(\d+)/;
		$start2 = $1;
		$end2 = $2;

		# no estamos imponiendo ningun valor de longitud
		# ¿que pasa con aquellos alineamiento sobrelapados?

		# usamos un valor arbitrario de 10 nt para ampliar la ventana de selección del
		# que será el mejor alineamiento
		if((($start1 - 10) < $start2 && $end1 > $end2)
		   || ($start1 < $start2 && ($end1 + 10) > $end2)){
		    # conservamos el segmento grande y eliminamos este..
		    delete($helitronMap{$key}{$key3});
		    goto REGRESO;
		} elsif((($start2 - 10) < $start1 && $end2 > $end1)
			|| ($start2 < $start1 && ($end2 + 10) > $end1)){
		    # tenemos un elemento más grande y tenemos que eliminar al anterior
		    delete($helitronMap{$key}{$key2});
		    goto REGRESO;
		}
	    }
	}
    }
}

####################################################################################################
# RECONSTRUCTING DATA ELEMENTS

for $key(sort keys %helitronMap){
    for $key2(sort keys %{$helitronMap{$key}}){
	$val = $helitronMapGhost{$key}{$key2};
	$helitronMap{$key}{$key2} =~ /(\w+)\|(\d+)\|(\d+\-\d+)\|(\d+\-\d+)/;
	$query    = $1;
	$tnumbers = $4;
	$tnumbers =~ /(\d+)\-(\d+)/;
	$startt   = $1;
	$endt     = $2;

	$segment{$query}{$key}{"id"}++;
	$id = $segment{$query}{$key}{"id"};

	$segment{$query}{$key}{$id} = $startt."-".$endt;  # guardamos las posiciones de inicio-fin
	# hash{target}{helitron}{index} = inicio.-.fin
	$backUPsegment{$query}{$key}{$id} = $val;  # mantenemos guardados los resultados del blast
	#backUPsegment{helitron}{indice}=genmapeo.|.longitud.|.inicio-finhelitron.|.inicio-fingenmapeo
    }
}

undef %helitronMap;
undef %helitronMapGhost;

####################################################################################################
## COMENZAMOS CON LA IDENTIFICACIÓN DE LOS GRUPOS DE CONSERVACIÓN

# corremos primero por cada elemento génico con el que se tiene homología... "key"
for $key(sort keys %segment){
    $index = 1;  # indicador del número de clases identificadas para cada gen "madre".
    # cada uno de los transposones que tiene un pedazo con este gen madre...
    for $TE(sort keys %{$segment{$key}}){   # comenzamos con cada transposón...
	# comenzamos a recorrer cada uno de los segmentos identificados para este gen transposon
	# $segment{$key}{$TE}{"id"}. contabiliza si hay más de un segmento en el mismo transposon que
	# captura un pedazo del mismo gen
	for($i=1; $i <= $segment{$key}{$TE}{"id"}; $i++){
	    $segment{$key}{$TE}{$i} =~ /(\d+)\-(\d+)/;
	    $start1 = $1;
	    $end1   = $2;
	    $val    = 0;   ## si tenemos dos variables que nos digan el valor de salida...
	    $refE   = 0;
	    for $key2(sort keys %{$Selement{$key}}){
		$Selement{$key}{$key2} =~ /(\d+||\d+\.\d+)\-(\d+||\d+\.\d+)/;
		$start2 = $1;
		$end2   = $2;
		# aquí comienza la identificación de lo segmentos...
		$val = &search_elem($start1,$end1,$start2,$end2);  # call the subroutine start, end
		# llaves de listElement
		$ref = 0;
		#if($val){
		# preguntamos si los elementos que el putativo nuevo segmento cumple las
		# características para estar dentro del grupo
		# $ref = &control_param($start1,$end1,$key,$key2);
		#}
		if($val){
		    $refE = 1;  # valor de cambio univalente indica que si el TE encontró un partner
		}

		next if(!$val);
		# val indica el tipo de relación encontrado, si el fragmento capturado por el
		# helitron a evaluar es ligeramente mayor que lo establecido por el grupo entonces:
		# val=1 y nos quedamos con las nuevas coordenadas..
		# si val=2 conservamos las mismas coordenadas,
		# si val=3 cambiamos el fin por las del nuevo elemento que estamos evaluando
		# y si val=4, entonces cambiamos el inicio

		# definimos promedios para las coordenadsa de inicio y fin...
		$startm = ($start1+$start2)/2;
		$endm   = ($end1+$end2)/2;

		if($val == 1){
		    $Selement{$key}{$key2}         = $startm."-".$endm;
		    $listElement{$key}{$key2}     .= "-".$TE.".".$i;
		    $backUPsegment{$key}{$TE}{$i} .= "|".$key2;
		} elsif($val == 2){
		    $Selement{$key}{$key2}         = $startm."-".$endm;
		    $listElement{$key}{$key2}     .= "-".$TE.".".$i;
		    $backUPsegment{$key}{$TE}{$i} .= "|".$key2;
		} elsif($val == 3){
		    $Selement{$key}{$key2}         = $startm."-".$endm;
		    $listElement{$key}{$key2}     .= "-".$TE.".".$i;
		    $backUPsegment{$key}{$TE}{$i} .= "|".$key2;
		} elsif($val == 4){
		    $Selement{$key}{$key2}         = $startm."-".$endm;
		    $listElement{$key}{$key2}     .= "-".$TE.".".$i;
		    $backUPsegment{$key}{$TE}{$i} .= "|".$key2;
		}
	    }
	    # definimos un nuevo grupo...
	    if(!$refE){
		$listElement{$key}{$index}     = $TE.".".$i;
		$Selement{$key}{$index}        = $start1."-".$end1;
		$backUPsegment{$key}{$TE}{$i} .= "|".$index;
		$index++;
	    }
	}
    }
}
## en teoría hasta aquí podemos definir/identificar los grupos de conservación.


# abrimos un archivo de salida para imprimir los datos
open(OUT, ">$ARGV[1]") or die "error i can not create the outfile\n";

print OUT "Gene\tClass\tHelitrons\tCapRange\tNoHelitrons\n";
for $key(sort keys %listElement){
    for $key2(sort {$a<=>$b} keys %{$listElement{$key}}){
	@data=split(/\-/,$listElement{$key}{$key2});
	#print OUT $key."\tclass\t".$key2."\t".$listElement{$key}{$key2}."\t".$Selement{$key}{$key2}."\t".scalar(@data)."\n";
	print OUT $key."\t".$key2."\t".$listElement{$key}{$key2}."\t".$Selement{$key}{$key2}."\t".scalar(@data)."\n";
    }
}
close(OUT);

open(OUT2, ">$ARGV[2]") or die "error i can not create the second outfile\n";
for $key(sort keys %backUPsegment){
    for $key2(sort keys %{$backUPsegment{$key}}){
	for $key3(sort keys %{$backUPsegment{$key}{$key2}}){
	    print OUT2 $key."\t".$key2."\t".$key3."\t".$backUPsegment{$key}{$key2}{$key3}."\n";
	}
    }
}
close(OUT2);
