import { LOGIC } from "../logic";
import { api } from "../sdkgen/api-generated";

api.fn.myBoards = (ctx) => {
	return LOGIC.boards.getUserBoards(ctx.request.extra.userId);
};

api.fn.getBoard = (ctx, { boardId }) => {
	return LOGIC.boards.getBoard(ctx.request.extra.userId, boardId);
};

api.fn.createBoard = async (ctx, { name, type }) => {
	await LOGIC.boards.createBoard(ctx.request.extra.userId, name, type);
};

api.fn.updateBoard = async (ctx, { boardId, name }) => {
	await LOGIC.boards.updateBoard(ctx.request.extra.userId, boardId, name);
};

api.fn.shareBoard = async (ctx, { email, boardId }) => {
	await LOGIC.boards.shareBoard(ctx.request.extra.userId, boardId, email);
};

api.fn.archiveBoard = async (ctx, { boardId }) => {
	await LOGIC.boards.archiveBoard(ctx.request.extra.userId, boardId);
};
