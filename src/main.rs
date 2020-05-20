use std::collections::HashMap;

pub struct SuffixTrie<'a>{
    pub text: String,
    pub suffix_trie: HashMap<char, &'a str>,
}

impl SuffixTrie<'_> {
    fn new(s: &str) -> Self {
        s.

    }

}

/*
 append '$' to t
set a hashmap to store a char --> remainder slice mapping
for i in length of t
curent dictionary = the hashmap
for character in t[i]
if character not in dictionary
add an empty vec
otherwise
set current_dictionary = current_dictionary[c]
*/

fn main() {
    println!("Hello, world!");
}
