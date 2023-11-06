# Use the official Julia image from Docker Hub
FROM julia:latest

# Set the working directory in the container
WORKDIR /env

# Copy only the project file(s) to the working directory
COPY Project.toml .
COPY Manifest.toml .

# Install project dependencies
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'

# Copy the entire project to the working directory
COPY . .

# Default command when the container starts
CMD ["julia", "--project=."]
