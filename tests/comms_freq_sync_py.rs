use comms_freq_sync_py::*;

#[test]
fn test_bit_reverse() {
    let correct_seqs: Vec<Vec<usize>> = vec![
        vec![0],
        vec![0,1],
        vec![0,2,1,3],
        vec![0,4,2,6,1,5,3,7],
        vec![0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15]
    ];
    for n_bits in 0..5 {
        println!("n_bits: {:?}", n_bits);
        for i in 0..n_bits {
            let rev = bit_reverse(i, n_bits);
            println!("rev: {:?} - corr: {:?}", rev, correct_seqs[n_bits as usize][i]);
            assert_eq!(rev, correct_seqs[n_bits as usize][i]);
        }
    }
}

// #[test]
