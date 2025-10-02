def detect_alphabet(S):
    return set(S)

def relative_frequency(S):
    freq = {}
    n = len(S)
    for char in S:
        freq[char] = freq.get(char, 0) + 1
    for char in freq:
        freq[char]/=n
    return freq

S = input('Input a sequence of 20 letters: ')

print('The alphabet of the sequence is: ', detect_alphabet(S))
print('The relative freq is: ', relative_frequency(S))