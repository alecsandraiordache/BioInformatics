import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def get_alignment_score(n1, n2, match=1, mismatch=-1, gap=0):
    if n1 == n2:
        return match
    elif n1 == '-' or n2 == '-':
        return gap
    else:
        return mismatch

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=0):
    n = len(seq1)
    m = len(seq2)

    score_matrix = np.zeros((n + 1, m + 1))
    
    for i in range(n + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(m + 1):
        score_matrix[0][j] = j * gap_penalty
        
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            
            score_matrix[i][j] = max(match, delete, insert)

    align1 = ""
    align2 = ""
    i, j = n, m
    
    path_x = [j]
    path_y = [i]

    while i > 0 and j > 0:
        score_current = score_matrix[i][j]
        score_diagonal = score_matrix[i - 1][j - 1]
        score_up = score_matrix[i - 1][j]
        score_left = score_matrix[i][j - 1]
        
        if score_current == score_diagonal + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_left + gap_penalty:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1
            
        path_x.append(j)
        path_y.append(i)

    while i > 0:
        align1 += seq1[i - 1]
        align2 += '-'
        i -= 1
        path_x.append(j)
        path_y.append(i)
    while j > 0:
        align1 += '-'
        align2 += seq2[j - 1]
        j -= 1
        path_x.append(j)
        path_y.append(i)

    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2, score_matrix, (path_x, path_y)

def calculate_metrics(align1, align2):
    matches = 0
    length = len(align1)
    match_line = ""
    
    for a, b in zip(align1, align2):
        if a == b and a != '-':
            matches += 1
            match_line += "|"
        elif a == '-' or b == '-':
            match_line += " "
        else:
            match_line += "." 
            
    similarity = (matches / length) * 100
    return matches, length, similarity, match_line

def visualize_alignment(score_matrix, path_coords, s1_label, s2_label):
    path_x, path_y = path_coords
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    sns.heatmap(score_matrix, cmap='magma', ax=ax1, cbar=True)
    ax1.set_title("Alignment Matrix Heatmap")
    ax1.set_xlabel("Sequence 2")
    ax1.set_ylabel("Sequence 1")
    
    rows, cols = score_matrix.shape
    path_grid = np.zeros((rows, cols))
    
    for x, y in zip(path_x, path_y):
        path_grid[y, x] = 1 
        
    from matplotlib.colors import ListedColormap
    cmap_custom = ListedColormap(['#FFFFE0', '#D32F2F']) 
  
    sns.heatmap(path_grid, cmap=cmap_custom, ax=ax2, linecolor='black', linewidths=0.5, cbar=False)
    ax2.set_title("Traceback Path")
    ax2.invert_yaxis() 
    
    ax2.invert_yaxis() 
    
    plt.tight_layout()
    plt.show()


S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"

GAP = 0
MATCH = 1
MISMATCH = -1

aligned_s1, aligned_s2, matrix, path = needleman_wunsch(S1, S2, MATCH, MISMATCH, GAP)
matches, length, similarity, visual_match = calculate_metrics(aligned_s1, aligned_s2)

print(f"S1: {S1}")
print(f"S2: {S2}")
print("Show Alignment:")
print(aligned_s1)
print(visual_match)
print(aligned_s2)
print(f"Matches    = {matches}")
print(f"Length     = {length}")
print(f"Similarity = {similarity:.0f} %")

visualize_alignment(matrix, path, S1, S2)