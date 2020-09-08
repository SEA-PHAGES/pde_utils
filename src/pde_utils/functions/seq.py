from pdm_utils.pipelines import export_db

from pde_utils.classes import kmers


def map_cds_to_seq(alchemist, cds_list, mol_type="dna", data_cache=None):
    if data_cache is None:
        data_cache = {}

    cds_to_seq = {}
    for cds in cds_list:
        if mol_type == "dna":
            parent_genome = data_cache.get(cds.genome_id)

            if parent_genome is None:
                parent_genome = export_db.get_single_genome(
                                            alchemist, cds.genome_id,
                                            data_cache=data_cache)

            cds.genome_length = parent_genome.length
            cds.set_seqfeature()
            cds.set_nucleotide_sequence(
                                parent_genome_seq=parent_genome.seq)

            cds_to_seq[cds] = str(cds.seq)
        elif mol_type == "amino":
            cds_to_seq[cds] = str(cds.translation)

    return cds_to_seq


def find_conserved_kmers(seq_id_to_sequence_map, kmer_length,
                         hash_count=None, fpp=0.0001):
    sequences = list(seq_id_to_sequence_map.values())
    num_sequences = len(sequences)
    approx_count = (len(sequences[0]) - kmer_length) * num_sequences

    cmsketch = kmers.CountMinSketch(approx_count,
                                    hash_count=hash_count, fpp=fpp)

    conserved_kmer_data = {}
    for seq_id, seq in seq_id_to_sequence_map.items():
        for i in range(len(seq) - kmer_length):
            subseq = seq[i:i+kmer_length]
            cmsketch.add(subseq)
            if cmsketch.check(subseq) >= num_sequences:
                kmer_data = conserved_kmer_data.get(subseq, [])
                kmer_data.append((seq_id, i))
                conserved_kmer_data[subseq] = kmer_data

    return conserved_kmer_data
