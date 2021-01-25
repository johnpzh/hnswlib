
data_base="/scratch/zpeng/data"
bin="./hnsw_search"
set -x

#### SIFT1M
data_dir="${data_base}/sift1m"
data="sift"
index="${data_dir}/${data}.hnsw"
query="${data_dir}/${data}_query.fvecs"
gt="${data_dir}/${data}.true-100_NN.v2.binary"
output="output.search_${data}.txt"
num_v=1
dim=128
num_q=10000
${bin} ${index} ${query} ${gt} ${num_v} ${dim} ${num_q} | tee ${output}

#### GIST1M
data_dir="${data_base}/gist1m"
data="gist"
index="${data_dir}/${data}.hnsw"
query="${data_dir}/${data}_query.fvecs"
gt="${data_dir}/${data}.true-100_NN.v2.binary"
output="output.search_${data}.txt"
num_v=1
dim=960
num_q=1000
${bin} ${index} ${query} ${gt} ${num_v} ${dim} ${num_q} | tee ${output}

#### DEEP10M
data_dir="${data_base}/deep1b"
data="deep10M"
index="${data_dir}/${data}.hnsw"
query="${data_dir}/${data}_query.fvecs"
gt="${data_dir}/${data}.true-100_NN.v2.binary"
output="output.search_${data}.txt"
num_v=10
dim=96
num_q=10000
${bin} ${index} ${query} ${gt} ${num_v} ${dim} ${num_q} | tee ${output}

#### SIFT100M
data_dir="${data_base}/sift1b"
data="sift100M"
index="${data_dir}/${data}.hnsw"
query="${data_dir}/${data}_query.fvecs"
gt="${data_dir}/${data}.true-100_NN.v2.binary"
output="output.search_${data}.txt"
num_v=100
dim=128
num_q=10000
${bin} ${index} ${query} ${gt} ${num_v} ${dim} ${num_q} | tee ${output}

#### DEEP100M
data_dir="${data_base}/deep1b"
data="deep100M"
index="${data_dir}/${data}.hnsw"
query="${data_dir}/${data}_query.fvecs"
gt="${data_dir}/${data}.true-100_NN.v2.binary"
output="output.search_${data}.txt"
num_v=100
dim=96
num_q=10000
${bin} ${index} ${query} ${gt} ${num_v} ${dim} ${num_q} | tee ${output}

set +x