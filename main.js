var grids = [];
var xn = 73;
var yn = 73;
var zn = 3;
var cell_count = xn*yn*zn;
var k = 0;
var total_density = 0;
var cur_total_density = 0;

function preload(){
}

function setup(){
	grids = grid.prototype.initializeGridArray(xn, yn, zn);
	grids[grid.prototype.getIdByChoords(12,12,1)].setRho(1.1);
	for (var i = 0; i < cell_count; i++) {
		total_density += grids[i].rho;
	}
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
	camera(0, 0, -300);
	background(51);
	translate(-1* xn* 5 ,-1*yn*5,0);

	fill(200,10,50);
	//sphere(100);
	cur_total_density = 0;
	for (var i = 0; i < cell_count; i+=1) {
		if(!grids[i].solid)
		grids[i].drawSphere(10, 10);

		cur_total_density +=grids[i].rho;
	};
	
	k+=1;
	console.log("draw",k);
	grid.prototype.runInstance(grids);//als callback??
	//console.log("total", total_density, "cur_total_density" , cur_total_density);
	console.log("density Change", (log(abs(total_density- cur_total_density))));
	if (total_density- cur_total_density < -5 ) noLoop();
}

