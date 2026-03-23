fn main() {
    println!("Hello, world!");
}

#[cfg(test)]
mod tests {
    use crate::main;
    #[test]
    fn check_main() {
        main()
    }
}
