BUILD_DIR ?= $(CURDIR)/../build/mbla-modern-lv1
PREFIX ?= $(BUILD_DIR)/install
BUILD_TYPE ?= Release
JOBS ?= 2

.PHONY: all configure test install

all: configure
	cmake --build "$(BUILD_DIR)" --config "$(BUILD_TYPE)" --parallel "$(JOBS)"

configure:
	cmake -S "$(CURDIR)" -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DBUILD_TESTING=ON

test: all
	cd "$(BUILD_DIR)" && ctest -C "$(BUILD_TYPE)" --output-on-failure

install: all
	cmake --install "$(BUILD_DIR)" --config "$(BUILD_TYPE)" --prefix "$(PREFIX)"
