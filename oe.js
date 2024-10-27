// Genetic Code Dictionary
const geneticCode = {
  'TTT': 'Phenylalanine', 'TTC': 'Phenylalanine', 'TTA': 'Leucine', 'TTG': 'Leucine',
  'CTT': 'Leucine', 'CTC': 'Leucine', 'CTA': 'Leucine', 'CTG': 'Leucine',
  'ATT': 'Isoleucine', 'ATC': 'Isoleucine', 'ATA': 'Isoleucine', 'ATG': 'Methionine',
  'GTT': 'Valine', 'GTC': 'Valine', 'GTA': 'Valine', 'GTG': 'Valine',
  'TCT': 'Serine', 'TCC': 'Serine', 'TCA': 'Serine', 'TCG': 'Serine',
  'CCT': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
  'ACT': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
  'GCT': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
  'TAT': 'Tyrosine', 'TAC': 'Tyrosine', 'TAA': 'Stop', 'TAG': 'Stop',
  'CAT': 'Histidine', 'CAC': 'Histidine', 'CAA': 'Glutamine', 'CAG': 'Glutamine',
  'AAT': 'Asparagine', 'AAC': 'Asparagine', 'AAA': 'Lysine', 'AAG': 'Lysine',
  'GAT': 'Aspartic Acid', 'GAC': 'Aspartic Acid', 'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
  'TGT': 'Cysteine', 'TGC': 'Cysteine', 'TGA': 'Stop', 'TGG': 'Tryptophan',
  'CGT': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
  'AGT': 'Serine', 'AGC': 'Serine', 'AGA': 'Arginine', 'AGG': 'Arginine',
  'GGT': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine'
};


// Function to translate a given DNA sequence to amino acids
function translateSequence() {
  const sequence = document.getElementById('sequence').value.toUpperCase();

  // Function to translate a sequence and stop at the stop codon
  function translate(sequence) {
    let aminoAcids = [];
    for (let i = 0; i < sequence.length; i += 3) {
      const codon = sequence.substring(i, i + 3);
      if (codon.length < 3) break; // Handle incomplete codon at the end
      if (geneticCode[codon]) {
        if (geneticCode[codon] === 'Stop') {
          aminoAcids.push('Stop');
          break; // Stop translation at the stop codon
        }
        aminoAcids.push(geneticCode[codon]);
      }
    }
    return aminoAcids;
  }

  // Translating the sequence
  const translatedAminoAcids = translate(sequence);

  // Displaying the amino acid sequence and stopping at the stop codon
  document.getElementById('output').innerText = `
    Amino Acid Sequence: ${translatedAminoAcids.join(' - ')}
  `;
}


// Function to count nucleotides
function countNucleotides() {
  const sequence = document.getElementById('sequence').value.toUpperCase();
  const counts = {
    'A': 0, 'T': 0, 'C': 0, 'G': 0, 'U': 0
  };

  for (let nucleotide of sequence) {
    if (counts.hasOwnProperty(nucleotide)) {
      counts[nucleotide]++;
    }
  }

  document.getElementById('output').innerText = `Nucleotide Counts: A = ${counts['A']}, T = ${counts['T']}, C = ${counts['C']}, G = ${counts['G']}, U = ${counts['U']}`;
}

// Function to find complementary strand
function findComplement() {
  const sequence = document.getElementById('sequence').value.toUpperCase();
  const complement = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'
  };

  let complementaryStrand = '';
  for (let nucleotide of sequence) {
    complementaryStrand += complement[nucleotide] || '';
  }

  document.getElementById('output').innerText = `Complementary Strand: ${complementaryStrand}`;
}

// Function to identify patterns
// function identifyPatterns() {
//   const sequence = document.getElementById('sequence').value.toUpperCase();
//   const palindromeRegex = /(?=((.)(.)(.?)\3\2))/g;
//   const palindromes = [];

//   let match;
//   while ((match = palindromeRegex.exec(sequence)) !== null) {
//     palindromes.push(match[0]);
//   }

//   document.getElementById('output').innerText = `Palindromic Patterns Found: ${palindromes.join(', ')}`;
// }

// Function to simulate mutations and identify mutation type
function simulateMutation() {
  const sequence = document.getElementById('sequence').value.toUpperCase();
  const mutationType = document.getElementById('mutation-type').value;
  const position = parseInt(document.getElementById('position').value, 10) - 1;
  const newNucleotide = document.getElementById('new-nucleotide').value.toUpperCase();

  let mutatedSequence = sequence.split('');

  // Applying the mutation based on the selected type
  if (mutationType === 'substitution') {
    mutatedSequence[position] = newNucleotide;
  } else if (mutationType === 'insertion') {
    mutatedSequence.splice(position, 0, newNucleotide);
  } else if (mutationType === 'deletion') {
    mutatedSequence.splice(position, 1);
  }

  // Joining the mutated sequence back into a string
  const mutatedSequenceStr = mutatedSequence.join('');

  // Function to translate a sequence and stop at the stop codon
  function translate(sequence) {
    let aminoAcids = [];
    for (let i = 0; i < sequence.length; i += 3) {
      const codon = sequence.substring(i, i + 3);
      if (geneticCode[codon]) {
        if (geneticCode[codon] === 'Stop') {
          aminoAcids.push('Stop');
          break; // Stop translation at the stop codon
        }
        aminoAcids.push(geneticCode[codon]);
      }
    }
    return aminoAcids;
  }

  // Translating the original and mutated sequences
  const originalAminoAcids = translate(sequence);
  const mutatedAminoAcids = translate(mutatedSequenceStr);

  // Determining the type of mutation
  let mutationDescription = "No mutation detected";
  if (originalAminoAcids.join(' - ') !== mutatedAminoAcids.join(' - ')) {
    if (mutatedAminoAcids.includes('Stop') && !originalAminoAcids.includes('Stop')) {
      mutationDescription = "Nonsense Mutation";
    } else if (originalAminoAcids.length === mutatedAminoAcids.length) {
      const differences = originalAminoAcids.filter((aa, index) => aa !== mutatedAminoAcids[index]);
      if (differences.length === 0) {
        mutationDescription = "Silent Mutation";
      } else {
        mutationDescription = "Missense Mutation";
      }
    } else {
      mutationDescription = "Frameshift Mutation";
    }
  } else {
    mutationDescription = "No significant change (Silent Mutation)";
  }

  // Displaying the mutated sequence, amino acid sequence, and mutation type
  document.getElementById('output').innerText = `
    Mutated Sequence: ${mutatedSequenceStr}
    Amino Acid Translation: ${mutatedAminoAcids.join(' - ') || 'Incomplete codon sequence'}
    Mutation Type: ${mutationDescription}
  `;
}



// Function to identify patterns
// function identifyPatterns() {
//   const sequence = document.getElementById('sequence').value.toUpperCase();

//   // Palindrome detection (existing logic)
//   const palindromeRegex = /(?=((.)(.)(.?)\3\2))/g;
//   const palindromes = [];
//   let match;
//   while ((match = palindromeRegex.exec(sequence)) !== null) {
//     palindromes.push(match[0]);
//   }


//   // Repeats detection: looking for sequences repeated twice consecutively
//   const repeats = [];
//   const repeatRegex = /([ATCGU]{2,})(?=\1)/g;
//   while ((match = repeatRegex.exec(sequence)) !== null) {
//     repeats.push(match[0]);
//   }

//   // Common motifs: Identifying specific biological motifs like "TATA" or "CGA"
//   const motifs = [];
//   const motifList = ['TATA', 'CGA', 'GATA', 'TTAA']; // Example motifs
//   motifList.forEach(motif => {
//     if (sequence.includes(motif)) {
//       motifs.push(motif);
//     }
//   });

//   // Consecutive nucleotide detection: finding sequences of four or more identical nucleotides
//   const consecutive = [];
//   const consecutiveRegex = /(A{4,}|T{4,}|C{4,}|G{4,}|U{4,})/g;
//   while ((match = consecutiveRegex.exec(sequence)) !== null) {
//     consecutive.push(match[0]);
//   }

//   // Displaying results
//   document.getElementById('output').innerText = `
//     Palindromic Patterns Found: ${palindromes.join(', ') || 'None'}
//     Repeated Patterns Found: ${repeats.join(', ') || 'None'}
//     Motifs Found: ${motifs.join(', ') || 'None'}
//     Consecutive Nucleotides: ${consecutive.join(', ') || 'None'}
//   `;
// }

function identifyPatterns() {
  const sequence = document.getElementById('sequence').value.toUpperCase();
  
  // Check for empty sequence
  if (!sequence) {
    document.getElementById('output').innerText = "Please enter a sequence.";
    return;
  }

  // Palindrome detection (improved logic)
  const minPalindromeLength = 4;
  const palindromes = [];
  for (let i = 0; i < sequence.length; i++) {
    for (let j = i + minPalindromeLength; j <= sequence.length; j++) {
      const subSequence = sequence.slice(i, j);
      if (isPalindrome(subSequence)) {
        palindromes.push(subSequence);
      }
    }
  }

  // Helper function to check if a sequence is a palindrome
  function isPalindrome(seq) {
    return seq === seq.split('').reverse().join('');
  }

  // Repeats detection: looking for sequences repeated twice consecutively
  const repeats = [];
  const repeatRegex = /([ATCGU]{2,})(?=\1)/g;
  let match;
  while ((match = repeatRegex.exec(sequence)) !== null) {
    repeats.push(match[0]);
  }

  // Common motifs: Identifying specific biological motifs like "TATA" or "CGA"
  const motifs = [];
  const motifList = ['TATA', 'CGA', 'GATA', 'TTAA']; // Example motifs
  motifList.forEach(motif => {
    if (sequence.includes(motif)) {
      motifs.push(motif);
    }
  });

  // Consecutive nucleotide detection: finding sequences of four or more identical nucleotides
  const consecutive = [];
  const consecutiveRegex = /(A{4,}|T{4,}|C{4,}|G{4,}|U{4,})/g;
  while ((match = consecutiveRegex.exec(sequence)) !== null) {
    consecutive.push(match[0]);
  }

  // Displaying results
  document.getElementById('output').innerText = `
    Palindromic Patterns Found: ${palindromes.join(', ') || 'None'}
    Repeated Patterns Found: ${repeats.join(', ') || 'None'}
    Motifs Found: ${motifs.join(', ') || 'None'}
    Consecutive Nucleotides: ${consecutive.join(', ') || 'None'}
  `;
}


