#!/usr/bin/env nextflow

//nextflow.enable.dsl=2
println "Prototype pipeline for analysis with Gentrius"


results_path="$PWD/results"

// Scripts and binaries ----------------------------------------------------------------------------------------
f_scripts="/Users/Olga/Projects/Science/Projects/terraces/gentrius_pipe_nf/scripts"
iqtree2_gentrius="$f_scripts/iqtree2_gentrius"
script_m_py="$f_scripts/script-gen_0_1_matrix.py"
script_sanity_check_2_r="$f_scripts/script-sanity_check_2-analyse_rf_all.r"
script_topology_summary="$f_scripts/"
// -------------------------------------------------------------------------------------------------------------
// Simulation parameters
taxa=[20]
genes=[5,10]
md=[30,50]
trees_num=2
print_lim=100
sanity_lim=10
// -------------------------------------------------------------------------------------------------------------
process random_matrix {

	input:
	each n from taxa
	each k from genes
	each m from md
	
	output:
	//stdout ch
	file "mrx_${n}_${k}_${m}_sub.*" into matrix

	script:
	"""
	file_m="mrx_${n}_${k}_${m}"
	python3 $script_m_py -n $n -k $k -m $m -o \$file_m
	t="0"
	while [ \$t -lt ${trees_num} ]
	do
		t=\$[\$t+1]
		cp \${file_m}_input \${file_m}_sub.\$t
	done
	"""	
}
// -------------------------------------------------------------------------------------------------------------
process random_tree {

	input:
	file mrx from matrix.flatten()

	output:
	tuple file("$mrx"),file("${mrx}.tree") into dataset_m_t

	script:
	tree_file="${mrx}.tree"
	"""
	n=`head -n 1 $mrx | awk -F " " '{print \$1}'`
	$iqtree2_gentrius -r \$n ${tree_file} -rlen 0 0 0 -quiet
        cat ${tree_file} | sed "s/:0.0000000000//g" | sed "s/T/sp/g" | sed "s/sp0/sp\${n}/g"> ${tree_file}-mod
        rm ${tree_file}
        mv ${tree_file}-mod ${tree_file}
        rm ${tree_file}.log
	"""

}


process gentrius_default_1MLN_10_MLN {

	input:
	tuple file(mrx),file(tree) from dataset_m_t

	output:
	tuple file("$mrx"),file("${tree}"),env(trees_num),env(stop_rule) into results_default
	file("${mrx}.res_default.log") into results_default_1MLN_10_MLN

	script:
	"""
	$iqtree2_gentrius -gentrius -pr_ab_matrix $mrx $tree -pre ${mrx}.res_default
	file="${mrx}.res_default.log"
	# ADD set ID to the log file
	echo "ADDED_INFO_SIMULATION_DATASET_ID: ${mrx.name}" >> \$file
	# Extract info about stopping rule
        stop_rule=0
        check_warning=`grep "stopping condition is active" \$file | wc -l | awk -F " " '{print \$1}'`
        if [ \$check_warning -eq 1 ]
        then
        check_type=`grep "Type of stopping rule: terrace size" \$file | wc -l | awk -F " " '{print \$1}'`
        	if [ \$check_type -eq 1 ]
        	then
        		stop_rule="1"
        	else
                	check_type=`cat \$file | grep "Type of stopping rule: number of visited intermediate trees" | wc -l | awk -F " " '{print \$1}'`
        		if [ \$check_type -eq 1 ]
                	then    
                		stop_rule="2"
                	else
                        	stop_rule="3"
                	fi
        	fi
	fi
	# ADD info about stopping rule to log
	echo "ADDED_INFO_STOP_RULE_ID: \${stop_rule}" >> \$file
	trees_num=`grep 'Number of trees on terrace' \$file | awk -F " " '{print \$6}'`

	"""
}


// ------------------------------------------------------------------------------------------------
// SPLIT CHANNEL WITH DEFAULT RESULTS FOR DOWNSTREAM ANALYSIS
// ------------------------------------------------------------------------------------------------
results_default
    .branch { 
	matrix,tree,tree_num, stop_rule ->
	x = Integer.valueOf(tree_num)
	y = Integer.valueOf(stop_rule)
	z = Integer.valueOf(print_lim)

        print_trees: x!=0 && x<z && y == 0
		return tuple(matrix,tree,tree_num)

	increase_t: x!=0 && (y==1 || y==2)
		return tuple(matrix,tree)

	complex_case: x==0 && y!=0
        	return tuple(matrix,tree)

	others: true
		return tuple(matrix,tree,x,y)
    } 
    .set { post_default } 
/*
post_default.print_trees.view{"$it print"}
post_default.increase_t.view{"$it increase"}
post_default.complex_case.view {"$it complex"}
post_default.others.view {"$it others"}
*/
// ------------------------------------------------------------------------------------------------

/*
# Split channel results_default into: with/without stop rule active
# Channel stop_rule=0 -> if trees_num in [1,100K] run with print -> in [2,100K] summarise trees
#								 -> in [1,5K]	run sanity check
# Channel stop_rule=1/2	-> trees_num == 0 -> run complex_analysis
#			-> trees_num != 0 -> run with mod thresholds: C1=100MLN C2=100MLN
# Channel stop_rule=3 -> nothing, ignore
*/


process summary_default {
	
	publishDir "$results_path/"

	input:
  	file("*.res_default.log") from results_default_1MLN_10_MLN.collect()
	
	output:
	//stdout ch_1
	file "results_all_summary_default" into summary_default

  	"""
	file_out="results_all_summary_default"
	if [ -e \$file_out ]
	then
		rm \$file_out
	fi
 	for file in *.res_default.log
	do
		# Extract info from ADDED
		datasetID=` grep 'ADDED_INFO_SIMULATION_DATASET_ID' \$file | awk -F " " '{print \$2}' `
		datasetID_short=`echo "\$datasetID" | awk -F "_sub" '{print \$1}'`
		stop_rule=` grep 'ADDED_INFO_STOP_RULE_ID' \$file | awk -F " " '{print \$2}' `
		# Extract info about input data
		taxon_num=` grep 'Number of taxa: ' \$file | awk -F " " '{print \$4}' `
		part_num=` grep 'Number of partitions' \$file | awk -F " " '{print \$4}' `
		md_percent=`grep 'missing entries in supermatrix' \$file | awk -F " " '{print \$7}'`
		uniq_taxon_num=` grep 'Number of special taxa' \$file | awk -F " " '{print \$9}' `
		taxon_num_init_tree=`grep 'Number of taxa on initial tree' \$file | awk -F " " '{print \$7}'`
		taxon_num_insert=`grep 'Number of taxa to be inserted' \$file | awk -F " " '{print \$7}'`

		# Extract info about results
		trees_num=`grep 'Number of trees on terrace' \$file | awk -F " " '{print \$6}'`
		trees_num_interm=`grep 'Number of intermediated trees visited' \$file | awk -F " " '{print \$6}'`
		dead_ends_num=`grep 'Number of dead ends encountered' \$file | awk -F " " '{print \$6}'`
		CPU=`grep 'Total CPU' \$file | awk -F " " '{print \$5, \$7}' `
		
		# Print results summary
		echo "\${file}__\$datasetID | TAXA \${taxon_num} PART \${part_num} MD \${md_percent} ROW_ZERO 0 COL_ZERO 1 UniqT \${uniq_taxon_num} ID \${datasetID_short} T_SIZE \${trees_num} CPU \${CPU} INT \${trees_num_interm} DEAD \${dead_ends_num} STOP_RULE \${stop_rule} | INIT_TREE \${taxon_num_init_tree} TAXA_TO_INSERT \${taxon_num_insert} MISS_PERCENT \${md_percent} | TERRAPHAST 0 0 0 0" >> \$file_out
	done
  	"""
}

process gentrius_print_trees {
c	input:
	tuple file(mrx),file(tree),val(tree_num) from post_default.print_trees

	output:
	tuple file(mrx),val(tree_num), file("$tree"), file("${mrx}.res_print.all_gen_terrace_trees") into generated_trees

	script:
	"""
	$iqtree2_gentrius -gentrius -pr_ab_matrix $mrx $tree -pre ${mrx}.res_print -g_print
	"""
}

// ------------------------------------------------------------------------------
// Split channel for downstream analysis: sanity and top_summary
// ------------------------------------------------------------------------------
generated_trees
    .branch {
	matrix,tree_num, tree, trees ->

	x=Integer.valueOf(tree_num)
	y=Integer.valueOf(sanity_lim)

	ch_sanity_check: x<y
		return tuple (matrix,tree,trees)
	ch_summary_trees: x>1
		return trees
	}

	.set{gen_t}

gen_t.ch_sanity_check.into{ch_sanity_1;ch_sanity_2;ch_sanity_3}
// ------------------------------------------------------------------------------

process sanity_check_1 {

	input:
	tuple file(mrx),file(tree),file(trees) from ch_sanity_1

	output:
	//tuple file(tree),env(check) into ch_sanity_check_results
	env(check) into ch_sanity_check_results

	script:
	"""
	# Check 1: an input tree must be in a set
	$iqtree2_gentrius -rf $tree $trees -pre ${tree}_init_tree -quiet
	file=${tree}_init_tree.rfdist
	check=`grep " 0" \$file | wc -l | awk -F " " '{print \$1}'`
	"""
}

process sanity_check_2 {
	input:
        tuple file(mrx),file(tree),file(trees) from ch_sanity_2
	
	output:
        //tuple file(tree),env(check) into ch_sanity_check_results_2
	env(check) into ch_sanity_check_results_2

	script:
        """
	# Check 2: all trees must be different
        $iqtree2_gentrius -rf_all $trees -pre ${trees}_all_distinct -quiet
	Rscript $script_sanity_check_2_r ${trees}_all_distinct.rfdist
	check=`grep "SUCCESS" ${trees}_all_distinct.rfdist.SANITY_RF_DIFF_RES | wc -l | awk -F " " '{print \$1}'`
	if [ "\$check" == "1" ]
	then
		rm ${trees}_all_distinct.rfdist.SANITY_RF_DIFF_RES
	fi
	"""
}

process sanity_check_3 {
        input:
        tuple file(mrx),file(tree),file(trees) from ch_sanity_3

        output:
        //tuple file(tree),env(check) into ch_sanity_check_results_3
	env(check) into ch_sanity_check_results_3

        script:
        """
        # Check 3: all trees must have the same set of subtrees
        $iqtree2_gentrius -gentrius -pr_ab_matrix $mrx $tree -g_query $trees -pre ${trees}_naive
	check=`grep "all trees belong to the same terrace" ${trees}_naive.log | wc -l | awk -F " " '{print \$1}'`
        """
}



//ch_sanity_check_results.mix(ch_sanity_check_results_2).mix(ch_sanity_check_results_3)

process sanity_check_results {

	publishDir "$results_path/"

	input:
	val(checks) from ch_sanity_check_results.mix(ch_sanity_check_results_2).mix(ch_sanity_check_results_3).collect()
	
	output:
	file("report_SANITY_CHECK") into results_sanity_check
	

	script:
	"""
	result=`cat checks* | grep "0" | wc -l | awk -F " " '{print \$1}'`
	if [ "\$result" -eq 0 ]
	then
		echo "SUCCESS: all datasets passed sanity checks" > report_SANITY_CHECK
	else
		echo "FAIL: not all datasets passed sanity checks!!!!" > report_SANITY_CHECK
	fi
	"""

}


process topological_summary {
	
	input:
	file(trees) from gen_t.ch_summary_trees

	output:
	file("${trees}_summary_topology") into summary_top_results
	
	script:
	"""
	# Build consensus tree
	\$iqtree2_gentrius -con $trees -minsup 0.999999 -quiet
	# Collect info:
	splits_remained=`grep "splits found" ${trees}.log | awk -F " " '{print \$1}'`
        splits_ignored=0
        check=`grep "discarded because frequency" ${trees}.log | wc -l | awk -F " " '{print \$1}'`
        if [ "\$check" -gt 0 ]
        then
        	splits_ignored=`grep "discarded because frequency" ${trees}.log | awk -F " " '{print \$1}'`
        fi
	trees_num=`wc -l $trees | awk -F " " '{print \$1}'`
	# Plot with multifurcating nodes and summary
	
	"""
}

