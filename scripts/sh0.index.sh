
data_base="/scratch/zpeng/data"
M=16

set -x

#### SIFT1M
data_dir="${data_base}/sift1m"
data="sift"
input="${data_dir}/${data}_base.fvecs"
output="${data_dir}/${data}.hnsw"
num_v=1
dim=128
efConstruction=200
./hnsw_index ${input} ${output} ${num_v} ${dim} ${M} ${efConstruction}

#### GIST1M
data_dir="${data_base}/gist1m"
data="gist"
input="${data_dir}/${data}_base.fvecs"
output="${data_dir}/${data}.hnsw"
num_v=1
dim=960
efConstruction=200
./hnsw_index ${input} ${output} ${num_v} ${dim} ${M} ${efConstruction}

#### DEEP10M
data_dir="${data_base}/deep1b"
data="deep10M"
input="${data_dir}/${data}_base.fvecs"
output="${data_dir}/${data}.hnsw"
num_v=10
dim=96
efConstruction=200
./hnsw_index ${input} ${output} ${num_v} ${dim} ${M} ${efConstruction}

#### SIFT100M
data_dir="${data_base}/sift1b"
data="sift100M"
input="${data_dir}/${data}_base.fvecs"
output="${data_dir}/${data}.hnsw"
num_v=100
dim=128
efConstruction=40
./hnsw_index ${input} ${output} ${num_v} ${dim} ${M} ${efConstruction}

#### DEEP100M
data_dir="${data_base}/deep1b"
data="deep100M"
input="${data_dir}/${data}_base.fvecs"
output="${data_dir}/${data}.hnsw"
num_v=100
dim=96
efConstruction=40
./hnsw_index ${input} ${output} ${num_v} ${dim} ${M} ${efConstruction}

set +x