all: build/CMakeCache.txt run_make
.PHONY: all

# Run cmake if not yet done or if CMakeLists.txt has changed.
build/CMakeCache.txt: CMakeLists.txt
	@echo "Running cmake"
	@mkdir -p build
	@cd build && cmake ..

run_make: build/CMakeCache.txt
	@echo "Running make"
	$(MAKE) -C build
.PHONY: run_make

update:
	@touch src/CMakeLists.txt
# 	@touch test/src/CMakeLists.txt
	$(MAKE) -C build
.PHONY: update

clean:
	@echo "Cleaning"
	@rm -rf build
# 	@rm -rf bin
# 	@rm -rf test/bin
.PHONY: clean
