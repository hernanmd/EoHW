Class {
	#name : #CGRehhCompos3r,
	#superclass : #CGRehhComposer,
	#category : #'EoHW-Core'
}

{ #category : #'instance creation' }
CGRehhCompos3r class >> exampleGenerate03 [

	self new
		outputRehhScript: 'rehh-comp.R';
		outputShellScript: 'run-rehh-comp-all-chr';
		outputRehhManhattanPlotScript: 'rehh-comp-manhattan-wgs.R';		
		cwd: 'D:\\Huellas_de_Seleccion\\Altiplano\\PHASED';
		chrRangeString: '1 29';
		firstCompPct: '5/8';
		secCompPct: '3/8';
		ancestralPopNr: 1;		
		
		popCount: 3;
		popNames: {  
			'Angus-chr' . 
			'Brangus-A-chr' };
		
		minPercGenoHap: 95;
		minPercGenoSnp: 50;
		ies2rsbMethod: 'bilateral';

		writeHistogramFunction: false;

		writeBifurcationDiagram: false;
		mrkFoc: 2429;
		allFoc: 2;
		nmrkL: 10;
		nmrkR: 10;
		bifurcationDiagramTitle: 'Bifurcation diagram (XX SNP on BTA1): Derived Allele';

		generate.
]

{ #category : #'code generation' }
CGRehhCompos3r >> generateMainRscript [
	" Private - Write files to disk "

	self outputRehhPatches asFileReference writeStreamDo: [ : stream |
		self writeRehhPatches: stream ].
	self outputRehhScript asFileReference writeStreamDo: [ : stream |
		self writeHeaderComment: stream.
		self writeLoadLibs: stream.
		self loadInputFiles: stream.
		self buildRehhFiles: stream.
		self writeRehhFunctions: stream.
		self writeHistogramFunctionTo: stream.
		self writeBifurcationDiagramFunction: stream.
		self writeRsbFunction: stream ].

]
