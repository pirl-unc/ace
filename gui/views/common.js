class Peptide {
  constructor(id, sequence) {
    this.id = id;
    this.sequence = sequence;
  }
}

class Assignment {
    constructor(peptide_id, peptide_sequence, pool_id, coverage_id, plate_id, well_id) {
        this.peptide_id = peptide_id;
        this.peptide_sequence = peptide_sequence;
        this.pool_id = pool_id;
        this.coverage_id = coverage_id;
        this.plate_id = plate_id;
        this.well_id = well_id;
    }
}

class AssignmentBenchReady {
    constructor(plate_id, well_id, peptide_ids, peptide_sequences) {
        this.plate_id = plate_id;
        this.well_id = well_id;
        this.peptide_ids = peptide_ids;
        this.peptide_sequences = peptide_sequences;
    }
}

class SpotCount {
    constructor(plate_id, well_id, spot_count) {
        this.plate_id = plate_id;
        this.well_id = well_id;
        this.spot_count = spot_count;
    }
}

class PreferredPeptidePair {
    constructor(peptide_1_id, peptide_1_sequence, peptide_2_id, peptide_2_sequence, similarity_score) {
        this.peptide_1_id = peptide_1_id;
        this.peptide_1_sequence = peptide_1_sequence;
        this.peptide_2_id = peptide_2_id;
        this.peptide_2_sequence = peptide_2_sequence;
        this.similarity_score = similarity_score;
    }
}

class DeconvolutionResult {
    constructor(peptide_id, peptide_sequence, hit_well_ids, hit_well_ids_count, peptide_spot_count, result) {
        this.peptide_id = peptide_id;
        this.peptide_sequence = peptide_sequence;
        this.hit_well_ids = hit_well_ids;
        this.hit_well_ids_count = hit_well_ids_count;
        this.peptide_spot_count = peptide_spot_count;
        this.result = result;
    }
}
