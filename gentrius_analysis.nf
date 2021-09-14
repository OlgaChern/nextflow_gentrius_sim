#!/usr/bin/env nextflow
// ------------------------------------------------------------------------------------------------
params.flag_plot=false
// ------------------------------------------------------------------------------------------------
println "Prototype pipeline for analysis with Gentrius"
// ------------------------------------------------------------------------------------------------
results_path="$PWD/results"
// ------------------------------------------------------------------------------------------------
// Scripts and binaries 
// ------------------------------------------------------------------------------------------------
iqtree2_gentrius="/Users/Olga/Projects/Science/Projects/terraces/nextflow_gentrius_sim/scripts/iqtree2_gentrius"
script_m_py="/Users/Olga/Projects/Science/Projects/terraces/nextflow_gentrius_sim/scripts/script-gen_0_1_matrix.py"
script_sanity_check_2_r="/Users/Olga/Projects/Science/Projects/terraces/nextflow_gentrius_sim/scripts/script-sanity_check_2-analyse_rf_all.r"
script_topology_summary="/Users/Olga/Projects/Science/Projects/terraces/nextflow_gentrius_sim/scripts/script-plot-contrees.r"
script_plot_matrix="/Users/Olga/Projects/Science/Projects/terraces/nextflow_gentrius_sim/scripts/script-plot-matrix.r"
// ------------------------------------------------------------------------------------------------
// PARAMETERS: SIMULATE DATASET
// ------------------------------------------------------------------------------------------------
// ---------------------------------------
// PARAMS MATRIX: major
// ---------------------------------------
taxa=[20]
genes=[10]
md=[30]
// ---------------------------------------
// PARAMS MATRIX: spread of missing data
// ---------------------------------------
r0_all=[1]	//#[0,1] # minimum number of 0 in each row
c0_all=[10]	//#[1,10,20,50]   # minimum number of 0 in each column
u_all=[25]	//#[0,5,20]    # % of rows with sum = 1


freq=[[1,10,60,30],[2,60,10,30]]		//[[100,0,0],[60,10,30],[10,60,30]]
probs=[[1,0.05,0.7,0.25],[2,0.1,0.6,0.3]]		//[[1,0,0],[0.1,0.6,0.3],[0.05,0.7,0.25]]

ch_freq=channel.from(freq)
ch_probs=channel.from(probs)
ch_freq.cross(ch_probs).into{r_info;c_info}


// ---------------------------------------
// Number of trees per matrix
// ---------------------------------------
trees_num=1
// ------------------------------------------------------------------------------------------------
// PARAMETERS: POST_ANALYSIS
// ------------------------------------------------------------------------------------------------
print_lim=100
sanity_lim=2
p_rm_leaves=0.10
// ---------------------------------------
// default stopping rule thresholds
// ---------------------------------------
stop_t_default=1000000		// C1: number of generated trees
stop_i_default=10000000		// C2: number of visited intermediate trees
// ---------------------------------------
// increased stopping rule thresholds
// ---------------------------------------
stop_t_increased=10000000	// C1: number of generated trees
stop_i_increased=10000000	// C2: number of visited intermediate trees
// ------------------------------------------------------------------------------------------------
// SIMULATE DATASET: RANDOM MATRIX of GENE PRESENCE_ABSENCE
// ------------------------------------------------------------------------------------------------
process random_matrix {

	publishDir "$results_path/matrices"

	input:
	each n from taxa
	each k from genes
	each m from md

	each r0 from r0_all
	each c0 from c0_all
	each u from u_all

	each r from r_info
        each c from c_info

	output:
	//stdout ch
	file "*sub.*" into matrix
	file "*_input" into plot_matrix

	script:
	id="${r[0][0]}${c[0][0]}"

        rf="${r[0][1]} ${r[0][2]} ${r[0][3]}"
        rp="${r[1][1]} ${r[1][2]} ${r[1][3]}"

        cf="${c[0][1]} ${c[0][2]} ${c[0][3]}"
        cp="${c[1][1]} ${c[1][2]} ${c[1][3]}"
	
	"""
	file_m="mrx_${n}_${k}_${m}_${r0}_${c0}_${u}_${id}"

	c0_num=\$[$n*$c0/100]
	u_num=\$[$n*$u/100]
	python3 $script_m_py -n $n -k $k -m $m -o \$file_m -rf $rf -rp $rp -cf $cf -cp $cp -r0 $r0 -c0 \$c0_num -u \$u_num
	t="0"
	while [ \$t -lt ${trees_num} ]
	do
		t=\$[\$t+1]
		cp \${file_m}_input \${file_m}_sub.\$t
	done

	"""	
}
//ch.view()
// ------------------------------------------------------------------------------------------------
// SIMULATE DATASET: RANDOM TREE
// ------------------------------------------------------------------------------------------------
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

// ------------------------------------------------------------------------------------------------
// ANALYSIS: MAIN GENTRIUS
// ------------------------------------------------------------------------------------------------
process gentrius {

	input:
	tuple file(mrx),file(tree) from dataset_m_t

	output:
	tuple file("$mrx"),file("${tree}"),env(trees_num),env(stop_rule) into results_gentrius_post
	file("${mrx}.res_gentrius.log") into results_gentrius_log

	script:
	"""
	$iqtree2_gentrius -gentrius -pr_ab_matrix $mrx $tree -pre ${mrx}.res_gentrius
	file="${mrx}.res_gentrius.log"
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
// SPLIT CHANNEL WITH GENTRIUS RESULTS FOR DOWNSTREAM ANALYSES
// ------------------------------------------------------------------------------------------------
// The resulting channels are exclusive
results_gentrius_post
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



//post_default.print_trees.view{"$it print"}
//post_default.increase_t.view{"$it increase"}
//post_default.complex_case.view {"$it complex"}
//post_default.others.view {"$it others"}

// ------------------------------------------------------------------------------------------------


//# Split channel results_default into: with/without stop rule active
//# Channel stop_rule=0 -> if trees_num in [1,100K] run with print -> in [2,100K] summarise trees
//#								 -> in [1,5K]	run sanity check
//# Channel stop_rule=1/2	-> trees_num == 0 -> run complex_analysis
//#			-> trees_num != 0 -> run with mod thresholds: C1=100MLN C2=100MLN
//# Channel stop_rule=3 -> nothing, ignore

// ------------------------------------------------------------------------------------------------
// POST_DEFAULT_ANALYSIS: increase threshold levels for stopping rules
// ------------------------------------------------------------------------------------------------
process gentrius_increased_t {

        input:
        tuple file(mrx),file(tree) from post_default.increase_t

        output:
        //tuple file("$mrx"),file("${tree}"),env(trees_num),env(stop_rule) into results_gentrius_post_t
        file("${mrx}.increased_t.res_gentrius.log") into results_gentrius_log_increased_t

        script:
        """
        $iqtree2_gentrius -gentrius -pr_ab_matrix $mrx $tree -pre ${mrx}.increased_t.res_gentrius -g_stop_t $stop_t_increased -g_stop_i $stop_i_increased
        file="${mrx}.increased_t.res_gentrius.log"
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



/*
// ------------------------------------------------------------------------------------------------
// SUMMARY: RESULTS GENTRIUS
// ------------------------------------------------------------------------------------------------

process summary_gentrius_log {
	
	publishDir "$results_path/"

	input:
  	file(log_files) from results_gentrius_log.concat(results_gentrius_log_increased_t,results_complex).collect()
	
	output:
	//stdout ch_1
	file "results_all_summary_default" into summary_default

  	"""
	file_out="results_all_summary_default"
	if [ -e \$file_out ]
	then
		rm \$file_out
	fi
 	for file in $log_files
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
*/

// ------------------------------------------------------------------------------------------------
// POST_DEFAULT_ANALYSIS: PRINT OUT GENERATED TREES
// ------------------------------------------------------------------------------------------------
process gentrius_print_trees {
	input:
	tuple file(mrx),file(tree),val(tree_num) from post_default.print_trees

	output:
	tuple file(mrx),val(tree_num), file("$tree"), file("${mrx}.res_print.all_gen_terrace_trees") into generated_trees

	script:
	"""
	$iqtree2_gentrius -gentrius -pr_ab_matrix $mrx $tree -pre ${mrx}.res_print -g_print
	"""
}

// ------------------------------------------------------------------------------
// COPY AND SPLIT channel for downstream analysis: sanity and top_summary
// ------------------------------------------------------------------------------
generated_trees.into {ch_sanity; ch_topology}

ch_sanity
    .branch {
	matrix,tree_num, tree, trees ->

	x=Integer.valueOf(tree_num)
	y=Integer.valueOf(sanity_lim)

	ch_sanity_check: x < y
		return tuple (matrix,tree,trees)
	}
	.set {gen_t}
	
gen_t.ch_sanity_check.into{ch_sanity_1;ch_sanity_2;ch_sanity_3}

ch_topology 
	.branch {
	matrix,tree_num, tree, trees ->

        x=Integer.valueOf(tree_num)
	ch_summary_trees: x > 1
		return trees
	}
	.set{ch_top}


// ------------------------------------------------------------------------------------------------
// POST_DEFAULT_ANALYSIS: SANITY CHECKs
// ------------------------------------------------------------------------------------------------
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


// ------------------------------------------------------------------------------------------------
// SUMMARY: SANITY CHECKs
// ------------------------------------------------------------------------------------------------
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
// ------------------------------------------------------------------------------------------------
// SUMMARY: TOPOLOGY
// ------------------------------------------------------------------------------------------------
process topological_summary {

	publishDir "$results_path/summary_topology/"
	
	input:
	file(trees) from ch_top.ch_summary_trees

	output:
	//file("${trees}_summary_topology") into summary_top_results
	file("${trees}.contree_plot_consensus_tree.pdf") into plots_consensus_trees
	
	script:
	"""
	# Build consensus tree
	$iqtree2_gentrius -con $trees -minsup 0.999999 -quiet
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
	Rscript $script_topology_summary ${trees}.contree \$splits_remained \$splits_ignored	
	"""
}
// ------------------------------------------------------------------------------------------------
// POST_DEFAULT_ANALYSIS: use REMOVE_LEAVES analysis to analyse complex datasets
// ------------------------------------------------------------------------------------------------
process analysis_complex_datasets {

	input:
	tuple file(mrx), file(tree) from post_default.complex_case

	output:
	file("${mrx}.rm_leaves.res_gentrius.log") into results_complex

	script:
	"""
	n=`head -n 1 $mrx | awk -F " " '{print \$1}'`
	rm_NUM=`echo "\$n*$p_rm_leaves" | bc | awk -F "." '{print \$1}'`

	$iqtree2_gentrius -gentrius -pr_ab_matrix $mrx $tree -pre ${mrx}.rm_leaves.res_gentrius -quiet -g_rm_leaves \$rm_NUM

	file="${mrx}.rm_leaves.res_gentrius.log"
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


if(params.flag_plot){


// ------------------------------------------------------------------------------------------------
// PLOT: MATRIX
// ------------------------------------------------------------------------------------------------
process plot_matrix {
	
	publishDir "$results_path/plots/matrices"

	input:
	file (mrx) from plot_matrix
	
	output:
	file ("*.pdf") into plotted_m
	
	script:
	"""
	Rscript $script_plot_matrix $mrx $mrx
	"""

}


}

// ------------------------------------------------------------------------------------------------
// SUMMARY: RESULTS GENTRIUS
// ------------------------------------------------------------------------------------------------

process summary_gentrius_log {

        publishDir "$results_path/"

        input:
        file(log_files) from results_gentrius_log.concat(results_gentrius_log_increased_t,results_complex).collect()

        output:
        //stdout ch_1
        file "results_all_summary_default" into summary_default

        """
        file_out="results_all_summary_default"
        if [ -e \$file_out ]
        then
                rm \$file_out
        fi
        for file in $log_files
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


