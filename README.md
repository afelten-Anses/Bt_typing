# Bt_typing

*Bacillus thuringiensis* serovars aizaiwai and kurstaki detection


Author: Arnaud Felten

Affiliation: [Ploufragan-Plouzané-Niort Laboratory, viral genetics and biosafety Unit, French Agency for Food, Environmental and Occupational Health (ANSES), 22440 Ploufragan, France](https://www.anses.fr/en/content/ploufragan-plouzan%C3%A9-niort-laboratory)

You can find the latest version of the tool at [https://github.com/afelten-Anses/Bt_typing](https://github.com/afelten-Anses/Bt_typing)


## External dependencies

* python3 (tested with 3.7.9)
* blast (tested with 2.10.1) 

**Note :** to use another blast version, you might have to remake the marker database using the *makeblastdb* command. 


## Parameters

Parameters of the script are available with one of its 3 options :

	Bt_detect.py
	Bt_detect.py -h
	Bt_detect.py --help
	
### Bt_detect parameters list

* -i : genome assemblie(s) path(s) in fasta format (REQUIRED)
* -t : table_bt.txt file path (default:*table_bt.txt*)
* -min_id : minimum percent of blast identity (default:90)
* -min_cov : minimum percent of blast coverage (default:90)
* -db_dir : marker sequences path (default:*db/markers*)
* -T : Number of threads to use (default:4)

## Test

Donwload or clone this git repository, then :

	cd Bt_typing
	python Bt_detect.py -i test/Bc_ATCC_14579.fasta test/Bt_kurstaki_HD1.fasta test/Leapi01.fasta -db_dir db/markers -t table_bt.txt > result.txt
