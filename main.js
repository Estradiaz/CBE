var grids = [];
var xn = 10;
var yn = 10;
var zn = 10;
var cell_count = xn*yn*zn;
var k = 0;

function preload(){
}

function setup(){
	grids = grid.prototype.initializeGridArray(xn, yn, zn);
	//grids.push(new grid(0,0,0,0));
	console.log(grids);
	createCanvas(600,600,WEBGL);
	background(51);
	//translate(width/2,height/2, -200);
	rotateZ(2);
}	

function draw(){
	var z = map(mouseX,0,2000,0,width);	
	z -= 1000;
	var x = map(mouseY,0,2000,0,height);
	//x -= 1000;
	camera(x, 0, z);
	background(51);
	translate(0,0,-400);
	//rotateY(k);
	
	fill(200,10,50);
	//sphere(100);
	for (var i = 0; i < cell_count; i+=1) {
		grids[i].drawSphere(10, 10);
	};
	
	k+=1;
	console.log("draw",k);
	grid.prototype.runInstance(grids);//als callback??

}
