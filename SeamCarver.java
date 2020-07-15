/*
Leon Deng
6/2/2019

SeamCarver takes as input a picture file (.png or .jpg) that should be cropped as well as the crop amounts
in horizontal and vertical pixels.
By deleting seams of pixels one at a time and by prioritizing lower-energy pixels, the resulting picture
will retain more of its distinguishing features than it would if it were traditionally cropped in chunks.

Sample Input:
Test1.png 100 100
Sample Output:
500-by-333 image
Time elapsed: 3.900801483 seconds.
new image size is 400 columns by 233 rows

Output images are attached as screenshots.
*/

import java.lang.IndexOutOfBoundsException;
import java.lang.NullPointerException;
import java.lang.IllegalArgumentException;

public class SeamCarver{

	private Picture picture;
	private double[][] energy;
	private int w;
	private int h;

	/*
	The constructor stores the inputted picture as well as the Picture's dimensions.
	Calculates the energy of each pixel in the current Picture and stores it in a 2d array.
	*/
	public SeamCarver(Picture picture){
		this.picture = picture;
		w = picture.width();
		h = picture.height();
		energy = new double[w][h];
		for(int i = 0; i < w; i ++){
			for(int j = 0; j < h; j++){
				energy[i][j] = energy(i,j);
			}
		}
	}

	//Returns the current picture
	public Picture picture(){
		return picture;
	}

	//Returns the width of current picture
	public int width(){
		return w;
	}

	//Returns the height of current picture
	public int height(){
		return h;
	}

	/*
	Uses the dual-gradient function to find the energy of each pixel.
	Finds the squared distance between the pixel's adjacent RGB colors in both the x- and y-directions.
	Adds these distances together to find the energy of the pixel.
	The dual-gradient function wraps around. For example, a pixel on the 0th row would be considered
	adjacent to a pixel on the (h-1)th row. Same with columns.
	*/
	public double energy(int x, int y){
		if (x < 0 || x > w - 1 || y < 0 || y > h - 1) {
            throw new IndexOutOfBoundsException();
        }

        int xc1 = (x - 1 + w) % w; //Eliminates need for conditionals, just wraps around
        int yc1 = y;
        int xc2 = (x + 1 + w) % w;
        int yc2 = y;

        int xc3 = x;
        int yc3 = (y - 1 + h) % h;
        int xc4 = x;
        int yc4 = (y + 1 + h) % h;

        int xRed = Math.abs(picture.get(xc1,yc1).getRed() - picture.get(xc2,yc2).getRed());
        int xGreen = Math.abs(picture.get(xc1,yc1).getGreen() - picture.get(xc2,yc2).getGreen());
        int xBlue = Math.abs(picture.get(xc1,yc1).getBlue() - picture.get(xc2,yc2).getBlue());

        int yRed = Math.abs(picture.get(xc3,yc3).getRed() - picture.get(xc4,yc4).getRed());
        int yGreen = Math.abs(picture.get(xc3,yc3).getGreen() - picture.get(xc4,yc4).getGreen());
        int yBlue = Math.abs(picture.get(xc3,yc3).getBlue() - picture.get(xc4,yc4).getBlue());

        double xGradient = Math.pow(xRed,2) + Math.pow(xGreen,2) + Math.pow(xBlue,2);
        double yGradient = Math.pow(yRed,2) + Math.pow(yGreen,2) + Math.pow(yBlue,2);

        return xGradient + yGradient;
	}

	/*
	Finds a 1-pixel-tall horizontal seam that sums up to the lowest energy.
	Imagines the starting point at the 0th column. Creates a 2d array corresponding to the pixels,
	then loops through each pixel and finds the lowest energy that could get to that pixel and stores
	that energy in the array. Also stores a 2d array for each pixel containing the vertical index of 
	the pixel that the lowest energy path came from. It's like a LinkedList.
	When the shortest paths (sps) array is done calculating, the lowest energy in the (w-1)th column is found.
	Then, using the from[] array, the pixel is traced back to the 0th column and the entire path is recorded
	and returned as the horizontal seam.
	*/
	public int[] findHorizontalSeam(){
		double[][] sps = new double[w][h];
		int[][] from = new int[w][h];

		for(int i = 0; i < h; i++){
			sps[0][i] = energy[0][i];
		}

		for(int i = 1; i < w; i++){
			for(int j = 0; j < h; j++){
				if(j == 0){
					if(sps[i-1][j] < sps[i-1][j+1]){
						sps[i][j] = sps[i-1][j] + energy[i][j];
						from[i][j] = j;
					}
					else{
						sps[i][j] = sps[i-1][j+1] + energy[i][j];
						from[i][j] = j+1;
					}
				}
				else if(j == h-1){
					if(sps[i-1][j-1] < sps[i-1][j]){
						sps[i][j] = sps[i-1][j-1] + energy[i][j];
						from[i][j] = j-1;
					}
					else{
						sps[i][j] = sps[i-1][j] + energy[i][j];
						from[i][j] = j;
					}
				}
				else{
					if(sps[i-1][j-1] <= sps[i-1][j] && sps[i-1][j-1] <= sps[i-1][j+1]){
						sps[i][j] = sps[i-1][j-1] + energy[i][j];
						from[i][j] = j-1;
					}
					else if(sps[i-1][j] <= sps[i-1][j-1] && sps[i-1][j] <= sps[i-1][j+1]){
						sps[i][j] = sps[i-1][j] + energy[i][j];
						from[i][j] = j;
					}
					else{
						sps[i][j] = sps[i-1][j+1] + energy[i][j];
						from[i][j] = j+1;
					}
				}
			}
		}

		int[] seam = new int[w];

		int min = 0;
		for(int i = 1; i < h; i++){
			if(sps[w-1][min] > sps[w-1][i]){
				min = i;
			}
		}
		seam[seam.length-1] = min;
		for(int i = w-1; i > 0; i--){
			seam[i-1] = from[i][seam[i]];
		}
		return seam;
	}

	// sequence of indices for vertical seam
	/*
	Dependent on EdgeWeightedDigraph, DijkstraSP, WeightedEdge
	Create EdgeWeightedDigraph of image, size width x height + 2: (0,0) is vertex 1, imaginary top vertex and bottom vertex
	for each vertex, add WeighedEdge from it to each of its 2/3 children, weight is just the energy of child
	Run pathTo() in DijkstraSP on the EdgeWeightedDigraph to produce the iterable
	Iterate through stack, then put them in array backwards, return the array
	
	public int[] findVerticalSeam(){
		EdgeWeightedDigraph edw = new EdgeWeightedDigraph(w * h + 2);

		//Connects top of Digraph to first row of pixels, and the bottom row to the bottom of Digraph
		for(int i = 0; i < w; i++){
			edw.addEdge(new DirectedEdge(0, i+1, energy[i][0]));
			edw.addEdge(new DirectedEdge(w*(h-1)+i+1, edw.V()-1, 0));
		}

		//Connects each pixel to bottom left, bottom, bottom right. Doesn't wrap around.
		for(int i = 0; i < w; i++){
			for(int j = 0; j < h - 1; j++){
				if(i == 0){
					edw.addEdge(new DirectedEdge(i+j*w+1, i+(j+1)*w+1, energy[i][j+1]));
					edw.addEdge(new DirectedEdge(i+j*w+1, i+(j+1)*w+2, energy[i+1][j+1]));
				}
				else if(i == w-1){
					edw.addEdge(new DirectedEdge(i+j*w+1, i+(j+1)*w, energy[i-1][j+1]));
					edw.addEdge(new DirectedEdge(i+j*w+1, i+(j+1)*w+1, energy[i][j+1]));
				}
				else{
					edw.addEdge(new DirectedEdge(i+j*w+1, i+(j+1)*w, energy[i-1][j+1]));
					edw.addEdge(new DirectedEdge(i+j*w+1, i+(j+1)*w+1, energy[i][j+1]));
					edw.addEdge(new DirectedEdge(i+j*w+1, i+(j+1)*w+2, energy[i+1][j+1]));
				}
			}
		}

		AcyclicSP sp = new AcyclicSP(edw,0);
		Stack<DirectedEdge> path = (Stack<DirectedEdge>) sp.pathTo(edw.V()-1);
		int[] output = new int[h];

		path.pop(); //Removes DirectedEdge from top of Digraph
		for(int i = 0; i < output.length; i++){
			output[i] = (path.pop().from() - 1) % w;
		}

		return output;
	}*/

	/*
	Removes the given horizontal seam from the Picture.
	Makes the Picture's height 1 less. The Picture is updated, as is the energy[][] list and h and w.
	Energy is lastly recalculated for the pixels that were affected by the seam, so a 2-pixel-tall
	section. This step is crucial - if the energy isn't recalculated after each step, removing more
	than 100 pixels or so would make the resultant picture look pretty strange.
	*/
	public void removeHorizontalSeam(int[] seam){
		if(seam == null){throw new NullPointerException();}
		if(seam.length != w){throw new IllegalArgumentException();}
		if(w == 1 || h == 1){throw new IllegalArgumentException();}

		Picture newPic = new Picture(w,h-1);
		double[][] newEnergy = new double[w][h-1];

		boolean passedSeam;
		for(int i = 0; i < w; i++){
			passedSeam = false;
			for(int j = 0; j < h; j++){
				if(j == seam[i]){
					passedSeam = true;
				}
				else if(passedSeam){
					newPic.set(i,j-1,picture.get(i,j));
					newEnergy[i][j-1] = energy[i][j];
				}
				else{
					newPic.set(i,j,picture.get(i,j));
					newEnergy[i][j] = energy[i][j];
				}
			}
		}
		picture = newPic;
		h--;

		int y1;
		int y2;
		for(int i = 0; i < w; i++){
			y1 = (seam[i] - 1 + h) % h;
			y2 = seam[i] % h;
			newEnergy[i][y1] = energy(i,y1);
			newEnergy[i][y2] = energy(i,y2);
		}

		energy = newEnergy;
	}


	//The vertical analog of findHorizontalSeam().
	public int[] findVerticalSeam(){

		double[][] sps = new double[w][h];
		int[][] from = new int[w][h];

		for(int i = 0; i < w; i++){
			sps[i][0] = energy[i][0];
		}

		for(int j = 1; j < h; j++){
			for(int i = 0; i < w; i++){
				if(i == 0){
					if(sps[i][j-1] < sps[i+1][j-1]){
						sps[i][j] = sps[i][j-1] + energy[i][j];
						from[i][j] = i;
					}
					else{
						sps[i][j] = sps[i+1][j-1] + energy[i][j];
						from[i][j] = i+1;
					}
				}
				else if(i == w-1){
					if(sps[i][j-1] < sps[i-1][j-1]){
						sps[i][j] = sps[i][j-1] + energy[i][j];
						from[i][j] = i;
					}
					else{
						sps[i][j] = sps[i-1][j-1] + energy[i][j];
						from[i][j] = i-1;
					}
				}
				else{
					if(sps[i-1][j-1] <= sps[i][j-1] && sps[i-1][j-1] <= sps[i+1][j-1]){
						sps[i][j] = sps[i-1][j-1] + energy[i][j];
						from[i][j] = i-1;
					}
					else if(sps[i][j-1] <= sps[i-1][j-1] && sps[i][j-1] <= sps[i+1][j-1]){
						sps[i][j] = sps[i][j-1] + energy[i][j];
						from[i][j] = i;
					}
					else{
						sps[i][j] = sps[i+1][j-1] + energy[i][j];
						from[i][j] = i+1;
					}
				}
			}
		}

		int[] seam = new int[h];

		int min = 0;
		for(int i = 1; i < w; i++){
			if(sps[min][h-1] > sps[i][h-1]){
				min = i;
			}
		}
		seam[seam.length-1] = min;
		for(int i = h-1; i > 0; i--){
			seam[i-1] = from[seam[i]][i];
		}
		return seam;
	}

	//The vertical analog of removeVerticalSeam().
	public void removeVerticalSeam(int[] seam){
		if(seam == null){throw new NullPointerException();}
		if(seam.length != h){throw new IllegalArgumentException();}
		if(w == 1 || h == 1){throw new IllegalArgumentException();}

		Picture newPic = new Picture(w-1,h);
		double[][] newEnergy = new double[w-1][h];

		boolean passedSeam;
		for(int j = 0; j < h; j++){
			passedSeam = false;
			for(int i = 0; i < w; i++){
				if(i == seam[j]){
					passedSeam = true;
				}
				else if(passedSeam){
					newPic.set(i-1,j,picture.get(i,j));
					newEnergy[i-1][j] = energy[i][j];
				}
				else{
					newPic.set(i,j,picture.get(i,j));
					newEnergy[i][j] = energy[i][j];
				}
			}
		}
		picture = newPic;
		w--;

		int x1;
		int x2;
		for(int i = 0; i < h; i++){
			x1 = (seam[i] - 1 + w) % w;
			x2 = seam[i] % w;
			newEnergy[x1][i] = energy(x1,i);
			newEnergy[x2][i] = energy(x2,i);
		}

		energy = newEnergy;
	}

	/*
	private void transpose(){
		Picture tPic = new Picture(h,w);
		double[][] tEnergy = new double[h][w];
		for(int i = 0; i < w; i++){
			for(int j = 0; j < h; j++){
				tPic.set(j,i,picture.get(i,j));
				tEnergy[j][i] = energy[i][j];
			}
		}
		energy = tEnergy;
		picture = tPic;
		int temp = w;
		w = h;
		h = temp;
	}*/

	public static void main(String[] args){
		if (args.length != 3) {
            StdOut.println("Usage:\njava ResizeDemo [image filename] [num columns to remove] [num rows to remove]");
            return;
        }

        Picture picture = new Picture(args[0]);
        int removeColumns = Integer.parseInt(args[1]);
        int removeRows = Integer.parseInt(args[2]); 

        StdOut.printf("%d-by-%d image\n", picture.width(), picture.height());
        SeamCarver sc = new SeamCarver(picture);

        long startTime = System.nanoTime();

        for (int i = 0; i < removeRows; i++) {
            int[] horizontalSeam = sc.findHorizontalSeam();
            sc.removeHorizontalSeam(horizontalSeam);
        }

        for (int i = 0; i < removeColumns; i++) {
            int[] verticalSeam = sc.findVerticalSeam();
            sc.removeVerticalSeam(verticalSeam);
        }

        long timeElapsed = System.nanoTime() - startTime;
        long conversion = 1000000000;

        System.out.println("Time elapsed: " + (double)timeElapsed/conversion + " seconds.");

        StdOut.printf("new image size is %d columns by %d rows\n", sc.width(), sc.height());
        picture.show();
        sc.picture().show();  
	}
}