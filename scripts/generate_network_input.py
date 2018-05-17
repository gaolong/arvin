import sys


snp_node_weights_file = sys.argv[1]
snp_gene_edge_weights_file = sys.argv[2]
gene_node_weights_file = sys.argv[3]
backbone_network_file = sys.argv[4]
output_edge_file = sys.argv[5]
output_node_file = sys.argv[6]


dp = {}
#f = file("snp_node_weights.disruption_log_p.txt")
f = file(snp_node_weights_file)
while True:
	line = f.readline()
	if len(line) == 0:
		break
	line = line.strip()
	mid = line.split("\t")
	dp[mid[0]] = mid[1]
f.close

hit = []
f = file(snp_gene_edge_weights_file)
fe = file(output_edge_file, "w")
fn = file(output_node_file, "w")
f.readline()
while True:
	line = f.readline()
	if len(line) == 0:
		break
	line = line.strip()
	mid = line.split("\t")	
	if mid[0] not in hit and mid[0] in dp:
		fe.write(line + "\tEP\n")
		fn.write(mid[0] + "\t" + dp[mid[0]] + "\teSNP\n")
		hit.append(mid[0])
f.close


f = file(gene_node_weights_file)
f.readline()
while True:
	line = f.readline()
	if len(line) == 0:
		break
	line = line.strip()
	mid = line.split("\t")
	if mid[0] not in hit:
		fn.write(line + "\tGene\n")
		hit.append(mid[0])
f.close
fn.close


f = file(backbone_network_file)
while True:
	line = f.readline()
	if len(line) == 0:
		break
	line = line.strip()
	mid = line.split("\t")
	if mid[-1] == "FI":
		fe.write(line + "\n")
f.close
fe.close

