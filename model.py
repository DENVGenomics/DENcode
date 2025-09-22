import math

def compute_transmission_probs(kij_matrix, gij_matrix, gamma=1.0):
    """
    Compute P(i→j) ∝ K_ij * exp(-gamma * G_ij), then normalize rows to get a probability matrix.
    """
    P_matrix = {}
    uuids_i = set(i for (i, j) in kij_matrix.keys())
    for i in uuids_i:
        row_scores = {}
        for j in uuids_i:
            if i == j:
                continue
            kij = kij_matrix.get((i, j), 0.0)
            gij = gij_matrix.get((i, j), 0.0)
            if kij > 0 and gij >= 0:
                row_scores[j] = kij * math.exp(-gamma * gij)
        total = sum(row_scores.values())
        if total > 0:
            P_matrix[i] = {j: score / total for j, score in row_scores.items()}
    return P_matrix
