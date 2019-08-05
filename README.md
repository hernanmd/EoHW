# Installation

```bash
 wget -O- get.pharo.org | bash
 ./pharo-ui Pharo.image &
```

From Pharo

```smalltalk
Metacello new
	baseline: 'BaselineOfEoHW';
	repository: 'github://hernanmd/BaselineOfEoHW/repository';
	onConflictUseIncoming;
	load 
```

# Examples

Open a browser on the main class:

```smalltalk
Smalltalk tools browser openOnClass: CGRehh class
```

and check the example methods on class side.


