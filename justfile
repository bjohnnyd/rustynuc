docker_cmd := "docker run --rm -it -v $(pwd):/drone/src -w /drone/src" 
# Command aliases
alias macos := build_apple-darwin
alias linux := build_unknown-linux-musl

build_apple-darwin:
	@echo "Compiling x86_64-apple-darwin target binary"
	mkdir -p target/binaries
	{{docker_cmd}} joseluisq/rust-linux-darwin-builder:latest \
	rustc -vV > target/binaries/latest.darwin.compilation.log
	{{docker_cmd}} \
	--env CC=o64-clang \
	--env CXX=o64-clang++ \
	--env LIBZ_SYS_STATIC=1 \
	joseluisq/rust-linux-darwin-builder:latest \
	cargo build --release --target x86_64-apple-darwin
	cd target/x86_64-apple-darwin/release && \
	tar cvzf ../../binaries/x86_64-apple-darwin.tar.gz ./rustynuc && \
	zip ../../binaries/x86_64-apple-darwin.zip ./rustynuc
	@echo "Finished compiling x86_64-apple-darwin target binary"

build_unknown-linux-musl:
	@echo "Compiling x86_64-unknown-linux-musl target binary"
	mkdir -p target/binaries
	{{docker_cmd}} joseluisq/rust-linux-darwin-builder:latest \
	rustc -vV > target/binaries/latest.linux.compilation.log
	{{docker_cmd}} joseluisq/rust-linux-darwin-builder:latest \
	cargo build --release --target x86_64-unknown-linux-musl
	cd target/x86_64-unknown-linux-musl/release  && \
	tar cvzf ../../binaries/x86_64-unknown-linux-musl.tar.gz ./rustynuc && \
	zip ../../binaries/x86_64-unknown-linux-musl.zip ./rustynuc
	@echo "Finished compiling x86_64-unknown-linux-musl target binary"

process:
	@echo "Running standard Testing, Linting and MSRV"
	cargo test -q --all
	cargo fmt -q --all -- --check
	cargo clippy --all-targets \
	--all-features -q -- \
	-D warnings \
	-A clippy::cognitive_complexity
