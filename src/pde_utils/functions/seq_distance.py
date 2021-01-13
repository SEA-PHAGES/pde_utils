def calculate_gcs(query_phams, target_phams):
    shared_phams = set(query_phams).intersection(set(target_phams))

    query_spp = len(shared_phams) / len(query_phams)
    target_spp = len(shared_phams) / len(target_phams)

    return (query_spp + target_spp) / 2
