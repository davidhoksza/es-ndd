FROM node:12

WORKDIR /usr/src/app

COPY package.json ./
COPY webpack.config.js ./
COPY .babelrc ./
COPY src ./src
COPY data/v2 ./data/v2

RUN npm install

RUN npm run build

CMD ["npx", "http-server", "./dist", "-p 8080"]

EXPOSE 8080



